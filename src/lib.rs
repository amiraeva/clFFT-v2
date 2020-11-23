// #[link(name = "clFFT", kind = "dylib")]
pub mod ffi;

fn init_setup_data() -> ffi::clfftSetupData {
    let mut sdata = ffi::clfftSetupData::default();

    sdata.major = ffi::clfftVersionMajor;
    sdata.minor = ffi::clfftVersionMinor;
    sdata.patch = ffi::clfftVersionPatch;
    sdata.debugFlags = 0;

    sdata
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clfft_trivial() -> ocl::Result<()> {
        const N: usize = 16;
        const DIM: ffi::clfftDim = ffi::clfftDim::CLFFT_1D;

        let mut err;

        let fft_setup = init_setup_data();

        unsafe {
            err = ffi::clfftSetup(&fft_setup);

            let mut phandle = ffi::clfftPlanHandle::default();
            let mut pro_que = ocl::ProQue::builder().build()?;
            let ctx = pro_que.context().as_core();
    
            err = ffi::clfftCreateDefaultPlan(&mut phandle, ctx.as_ptr(), DIM, &N);
            err = ffi::clfftSetPlanPrecision(phandle, ffi::clfftPrecision::CLFFT_SINGLE_FAST);
        }


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

        println!("The value at index [{}] is now '{}'!", 200007, read_buf[200007]);

        Ok(())
    }
}
