use super::*;

use ocl::{Buffer, MemFlags};

#[test]
fn clfft_trivial() -> Result<()> {
    // Prepare some data
    let mut source = vec![0.0; 128];
    for i in 0..source.len() / 2 {
        let x = std::f64::consts::PI * 4.0 * (i as f64 / 2.0 / source.len() as f64);
        source[2 * i] = x.sin();
        source[2 * i + 1] = x.sin();
    }

    let ctx = ocl::Context::builder().build()?;
    let device = ocl::Device::first(ocl::Platform::first()?)?;
    let queue = ocl::Queue::new(&ctx, device, None)?;

    // Create buffers
    let mut in_buffer = unsafe {
        Buffer::builder()
            .queue(queue.clone())
            .flags(MemFlags::new().read_write())
            .len(source.len())
            .use_host_slice(&source)
            .build()?
    };

    let mut res_buffer = Buffer::<f64>::builder()
        .queue(queue.clone())
        .flags(MemFlags::new().write_only())
        .len(2 * source.len())
        .build()
        .expect("Failed to create GPU result buffer");

    // Make a plan
    let mut plan = builder::<f64>()?
        .precision(Precision::Precise)
        .dims(source.len())
        .input_layout(Layout::Real)
        .output_layout(Layout::HermitianInterleaved)
        .bake_out_of_place_plan(&queue, &ctx)?;

    // Execute plan
    plan.enq(Direction::Forward, &mut in_buffer, &mut res_buffer)
        .unwrap();

    // Wait for calculation to finish and read results
    res_buffer
        .cmd()
        .read(&mut source)
        .enq()
        .expect("Transferring result vector from the GPU back to memory failed");

    queue.finish()?;
    println!("{:?}", source);

    Ok(())
}

#[test]
fn clfft_trivial2() -> Result<()> {
    const N: usize = 16;

    let ctx = ocl::Context::builder().build()?;
    let device = ocl::Device::first(ocl::Platform::first()?)?;
    let queue = ocl::Queue::new(&ctx, device, None)?;

    let x = vec![0f32; 2 * N];
    let mut bufx = unsafe {
        ocl::Buffer::builder()
            .queue(queue.clone())
            .flags(ocl::MemFlags::new().read_write())
            .len(2 * N)
            .use_host_slice(&x)
            .build()?
    };

    bufx.cmd().write(&x).enq()?;

    // // Make a plan
    let mut plan = builder::<f32>()?
        .precision(Precision::Precise)
        .dims(N)
        .input_layout(Layout::ComplexInterleaved)
        .output_layout(Layout::ComplexInterleaved)
        .bake_inplace_plan(&queue, &ctx)?;

    plan.enq(Direction::Forward, &mut bufx)?;

    Ok(())
}

const SRC: &str = r#"
    __kernel void add(__global float* buffer, float scalar) {
        buffer[get_global_id(0)] += scalar;
    }
"#;

#[test]
fn trivial() -> ocl::Result<()> {
    let proque = ocl::ProQue::builder().src(SRC).dims(1 << 20).build()?;
    let buf = proque.create_buffer::<f32>()?;

    let kernel = proque
        .kernel_builder("add")
        .arg(&buf)
        .arg(10.0f32)
        .build()?;

    unsafe { kernel.enq()? };

    let mut read_buf: Box<[f32]> = (0..buf.len()).map(|_| 0.0).collect();

    buf.read(&mut *read_buf).enq()?;

    println!(
        "The value at index [{}] is now '{}'!",
        200007, read_buf[200007]
    );

    Ok(())
}
