use std::marker::PhantomData;

pub mod ffi;

fn init_setup_data() -> ffi::clfftSetupData {
    let mut sdata = ffi::clfftSetupData::default();

    sdata.major = ffi::clfftVersionMajor;
    sdata.minor = ffi::clfftVersionMinor;
    sdata.patch = ffi::clfftVersionPatch;
    sdata.debugFlags = 0;

    sdata
}
#[derive(Debug)]
pub enum Error {
    Ffi(ffi::clfftStatus),
    Ocl(ocl::Error),
    UnspecifiedDimensionError,
}

impl From<ffi::clfftStatus> for Result<()> {
    fn from(val: ffi::clfftStatus) -> Self {
        match val {
            ffi::clfftStatus::CLFFT_SUCCESS => Ok(()),
            _ => Err(Error::Ffi(val)),
        }
    }
}

impl From<ocl::Error> for Error {
    fn from(val: ocl::Error) -> Self {
        Self::Ocl(val)
    }
}

type Result<T> = std::result::Result<T, Error>;

macro_rules! clfft_try {
    ( $result_expr: expr) => {
        let err = $result_expr;
        Result::from(err)?;
    }
}

macro_rules! clfft_panic {
    ( $result_expr: expr) => {
        let result = $result_expr;
        if result != ffi::clfftStatus::CLFFT_SUCCESS {
            panic!("Unexpeced error in foreign library: {:?}", result);
        }
    };
}

pub struct SetupData(ffi::clfftSetupData);

impl SetupData {
    fn new() -> Result<Self> {
        let sdata = init_setup_data();
        let err = unsafe { ffi::clfftSetup(&sdata) };
        Result::from(err)?;
        Ok(Self(sdata))
    }
}

impl Drop for SetupData {
    fn drop(&mut self) {
        let err = unsafe { ffi::clfftTeardown() };
        if let Err(e) = Result::from(err) {
            panic!("Error when trying to drop clFFT resources\n{:?}", e);
        }
    }
}

/// A trait for all paremeters supported by `clFFT`.
pub trait ClFftPrm: ocl::OclPrm {
    /// Is the type a double precision type.
    fn is_dbl_precision() -> bool;
}

impl ClFftPrm for f32 {
    fn is_dbl_precision() -> bool {
        false
    }
}

impl ClFftPrm for f64 {
    fn is_dbl_precision() -> bool {
        true
    }
}

fn translate_to_fft_dim(dims: ocl::SpatialDims) -> ffi::clfftDim {
    match dims.dim_count() {
        1 => ffi::clfftDim::CLFFT_1D,
        2 => ffi::clfftDim::CLFFT_2D,
        3 => ffi::clfftDim::CLFFT_3D,
        n => panic!("Number of dimensions must be 1, 2 or 3, but it is {}", n),
    }
}

/// Specify the expected precision of each FFT.
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Precision {
    Precise,
    Fast,
}

fn translate_precision<T: ClFftPrm>(precision: Precision) -> ffi::clfftPrecision {
    let is_f64 = T::is_dbl_precision();
    match precision {
        Precision::Precise => {
            if is_f64 {
                ffi::clfftPrecision::CLFFT_DOUBLE
            } else {
                ffi::clfftPrecision::CLFFT_SINGLE
            }
        }
        Precision::Fast => {
            if is_f64 {
                ffi::clfftPrecision::CLFFT_DOUBLE_FAST
            } else {
                ffi::clfftPrecision::CLFFT_SINGLE_FAST
            }
        }
    }
}

fn translate_precision_back(precision: ffi::clfftPrecision) -> Precision {
    match precision {
        ffi::clfftPrecision::CLFFT_SINGLE => Precision::Precise,
        ffi::clfftPrecision::CLFFT_DOUBLE => Precision::Precise,
        ffi::clfftPrecision::CLFFT_SINGLE_FAST => Precision::Fast,
        ffi::clfftPrecision::CLFFT_DOUBLE_FAST => Precision::Fast,
        ffi::clfftPrecision::ENDPRECISION => panic!("ENDPRECISION should never be returned"),
    }
}

/// Specify the expected layouts of the buffers.
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Layout {
    ComplexInterleaved,
    ComplexPlanar,
    HermitianInterleaved,
    HermitianPlanar,
    Real,
}

fn translate_layout(layout: Layout) -> ffi::clfftLayout {
    match layout {
        Layout::ComplexInterleaved => ffi::clfftLayout::CLFFT_COMPLEX_INTERLEAVED,
        Layout::ComplexPlanar => ffi::clfftLayout::CLFFT_COMPLEX_PLANAR,
        Layout::HermitianInterleaved => ffi::clfftLayout::CLFFT_HERMITIAN_INTERLEAVED,
        Layout::HermitianPlanar => ffi::clfftLayout::CLFFT_HERMITIAN_PLANAR,
        Layout::Real => ffi::clfftLayout::CLFFT_REAL,
    }
}

fn translate_layout_back(layout: ffi::clfftLayout) -> Layout {
    match layout {
        ffi::clfftLayout::CLFFT_COMPLEX_INTERLEAVED => Layout::ComplexInterleaved,
        ffi::clfftLayout::CLFFT_COMPLEX_PLANAR => Layout::ComplexPlanar,
        ffi::clfftLayout::CLFFT_HERMITIAN_INTERLEAVED => Layout::HermitianInterleaved,
        ffi::clfftLayout::CLFFT_HERMITIAN_PLANAR => Layout::HermitianPlanar,
        ffi::clfftLayout::CLFFT_REAL => Layout::Real,
        ffi::clfftLayout::ENDLAYOUT => panic!("ENDLAYOUT should never be returned"),
    }
}

/// Specify the expected direction of each FFT, time or the frequency domains
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Direction {
    Forward,
    Backward,
}

fn translate_direction(direction: Direction) -> ffi::clfftDirection {
    match direction {
        Direction::Forward => ffi::clfftDirection::CLFFT_FORWARD,
        Direction::Backward => ffi::clfftDirection::CLFFT_BACKWARD,
    }
}

/// pecify wheter the input buffers are overwritten with results
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Location {
    Inplace,
    OutOfPlace,
}

fn translate_location(location: Location) -> ffi::clfftResultLocation {
    match location {
        Location::Inplace => ffi::clfftResultLocation::CLFFT_INPLACE,
        Location::OutOfPlace => ffi::clfftResultLocation::CLFFT_OUTOFPLACE,
    }
}

fn translate_location_back(location: ffi::clfftResultLocation) -> Location {
    match location {
        ffi::clfftResultLocation::CLFFT_INPLACE => Location::Inplace,
        ffi::clfftResultLocation::CLFFT_OUTOFPLACE => Location::OutOfPlace,
        ffi::clfftResultLocation::ENDPLACE => panic!("ENDPLACE should never be returned"),
    }
}

/// Builder for a FFT plan.
pub struct FftPlanBuilder<'a, T: ClFftPrm> {
    data_type: PhantomData<T>,
    precision: Precision,
    dims: Option<ocl::SpatialDims>,
    input_layout: Layout,
    output_layout: Layout,
    forward_scale: Option<f32>,
    backward_scale: Option<f32>,
    batch_size: Option<usize>,
    setup_data: PhantomData<&'a SetupData>,
}

/// Creates a builder for baking a new FFT plan.
pub fn builder<T: ClFftPrm>(_: &SetupData) -> FftPlanBuilder<'_, T> {
    FftPlanBuilder {
        data_type: PhantomData,
        precision: Precision::Precise,
        dims: None,
        input_layout: Layout::ComplexInterleaved,
        output_layout: Layout::ComplexInterleaved,
        forward_scale: None,
        backward_scale: None,
        batch_size: None,
        setup_data: PhantomData,
    }
}

impl<T: ClFftPrm> FftPlanBuilder<'_, T> {
    /// Set the floating point precision of the FFT data.
    pub fn precision(mut self, precision: Precision) -> Self {
        self.precision = precision;
        self
    }

    /// Set the expected layout of the input buffer.
    pub fn input_layout(mut self, input_layout: Layout) -> Self {
        self.input_layout = input_layout;
        self
    }

    /// Set the expected layout of the output buffer.
    pub fn output_layout(mut self, output_layout: Layout) -> Self {
        self.output_layout = output_layout;
        self
    }

    /// Set the dimensionality of the FFT transform; describes how many elements are in the array.
    ///
    /// If the data is complex then the dimesions are specified are per complex number. In practice that
    /// means that the dimensions should be half the size of the buffers.
    pub fn dims<D: Into<ocl::SpatialDims>>(mut self, dims: D) -> Self {
        self.dims = Some(dims.into());
        self
    }

    /// Set the scaling factor that is applied to the FFT data.
    pub fn forward_scale(mut self, scale: f32) -> Self {
        self.forward_scale = Some(scale);
        self
    }

    /// Set the scaling factor that is applied to the FFT data.
    pub fn backward_scale(mut self, scale: f32) -> Self {
        self.backward_scale = Some(scale);
        self
    }

    /// Set the number of discrete arrays that the plan can concurrently handle.
    pub fn batch_size(mut self, scale: usize) -> Self {
        self.batch_size = Some(scale);
        self
    }

    /// Creates a plan for an inplace FFT.
    pub fn bake_inplace_plan<'a>(
        &mut self,
        queue: &'a ocl::Queue,
        ctx: &ocl::Context,
    ) -> Result<FftInplacePlan<'a, T>> {
        let handle = bake_plan::<T>(self, queue, ctx, Location::Inplace)?;
        Ok(FftInplacePlan {
            handle: handle,
            queue,
            data_type: std::marker::PhantomData,
            wait_list: None,
            dest_list: None,
        })
    }

    /// Creates a plan for an out of place FFT.
    pub fn bake_out_of_place_plan<'a>(
        &mut self,
        queue: &'a ocl::Queue,
        ctx: &ocl::Context,
    ) -> Result<FftOutOfPlacePlan<'a, T>> {
        let handle = bake_plan::<T>(self, queue, ctx, Location::OutOfPlace)?;
        Ok(FftOutOfPlacePlan {
            handle: handle,
            queue,
            data_type: std::marker::PhantomData,
            wait_list: None,
            dest_list: None,
        })
    }
}

fn bake_plan<T: ClFftPrm>(
    builder: &mut FftPlanBuilder<T>,
    queue: &ocl::Queue,
    ctx: &ocl::Context,
    location: Location,
) -> Result<ffi::clfftPlanHandle> {
    let context = ctx.as_ptr();
    let dims = builder.dims.unwrap();

    let dim = translate_to_fft_dim(dims);
    let lengths = dims
        .to_lens()
        .map_err(|_| Error::UnspecifiedDimensionError)?;

    let mut plan: ffi::clfftPlanHandle = 0;
    clfft_try!(unsafe { ffi::clfftCreateDefaultPlan(&mut plan, context, dim, lengths.as_ptr()) });
    let precision = translate_precision::<T>(builder.precision);
    clfft_try!(unsafe { ffi::clfftSetPlanPrecision(plan, precision) });
    let input_layout = translate_layout(builder.input_layout);
    let output_layout = translate_layout(builder.output_layout);
    clfft_try!(unsafe { ffi::clfftSetLayout(plan, input_layout, output_layout) });
    match builder.forward_scale {
        None => (),
        Some(s) => {
            clfft_try!(unsafe {
                ffi::clfftSetPlanScale(plan, ffi::clfftDirection::CLFFT_FORWARD, s)
            });
        }
    }
    match builder.backward_scale {
        None => (),
        Some(s) => {
            clfft_try!(unsafe {
                ffi::clfftSetPlanScale(plan, ffi::clfftDirection::CLFFT_BACKWARD, s)
            });
        }
    }
    let location = translate_location(location);
    clfft_try!(unsafe { ffi::clfftSetResultLocation(plan, location) });
    match builder.batch_size {
        None => (), // Use default
        Some(s) => {
            clfft_try!(unsafe { ffi::clfftSetPlanBatchSize(plan, s) });
        }
    }

    let mut queue = queue.as_ptr();
    clfft_try!(unsafe {
        ffi::clfftBakePlan(plan, 1, &mut queue, None, 0 as *mut ::std::os::raw::c_void)
    });
    Ok(plan)
}

/// A plan is a repository of state for calculating FFT's.  Allows the runtime to pre-calculate kernels, programs
/// and buffers and associate them with buffers of specified dimensions.
pub struct FftInplacePlan<'a, T: ClFftPrm> {
    handle: ffi::clfftPlanHandle,
    queue: &'a ocl::Queue,
    data_type: std::marker::PhantomData<T>,
    wait_list: Option<&'a ocl::EventList>,
    dest_list: Option<&'a mut ocl::EventList>,
}

/// A plan is a repository of state for calculating FFT's.  Allows the runtime to pre-calculate kernels, programs
/// and buffers and associate them with buffers of specified dimensions.
pub struct FftOutOfPlacePlan<'a, T: ClFftPrm> {
    handle: ffi::clfftPlanHandle,
    queue: &'a ocl::Queue,
    data_type: std::marker::PhantomData<T>,
    wait_list: Option<&'a ocl::EventList>,
    dest_list: Option<&'a mut ocl::EventList>,
}

impl<'a, T: ClFftPrm> FftInplacePlan<'a, T> {
    /// Specifies the list of events to wait on before the command will run.
    pub fn ewait(mut self, wait_list: &'a ocl::EventList) -> Self {
        self.wait_list = Some(wait_list);
        self
    }

    /// Specifies the destination list or empty event for a new, optionally
    /// created event associated with this command.
    pub fn enew(mut self, new_event_dest: &'a mut ocl::EventList) -> Self {
        self.dest_list = Some(new_event_dest);
        self
    }

    /// Enqueues the FFT so that it gets performed on the device.
    pub fn enq(&mut self, direction: Direction, buffer: &mut ocl::Buffer<T>) -> Result<()> {
        let input_len = self.dims().to_len()
            * if self.input_layout() == Layout::Real {
                1
            } else {
                2
            };

        if input_len != buffer.len() {
            // return Err(format!("FFT plan requires that input buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", input_len).into());
            // TODO: convert to error
            panic!("FFT plan requires that input buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", input_len);
        }

        enqueue::<T>(
            self.handle,
            direction,
            self.queue,
            buffer,
            None,
            &self.wait_list,
            &mut self.dest_list,
        )
    }
}

impl<'a, T: ClFftPrm> FftOutOfPlacePlan<'a, T> {
    /// Specifies the list of events to wait on before the command will run.
    pub fn ewait(mut self, wait_list: &'a ocl::EventList) -> Self {
        self.wait_list = Some(wait_list);
        self
    }

    /// Specifies the destination list or empty event for a new, optionally
    /// created event associated with this command.
    pub fn enew(mut self, new_event_dest: &'a mut ocl::EventList) -> Self {
        self.dest_list = Some(new_event_dest);
        self
    }

    /// Enqueues the FFT so that it gets performed on the device.
    pub fn enq(
        &mut self,
        direction: Direction,
        buffer: &ocl::Buffer<T>,
        result: &mut ocl::Buffer<T>,
    ) -> Result<()> {
        let input_len = self.dims().to_len()
            * if self.input_layout() == Layout::Real {
                1
            } else {
                2
            };

        if input_len != buffer.len() {
            // return Err(format!("FFT plan requires that input buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", input_len).into());
            // TODO: convert to error
            panic!("FFT plan requires that input buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", input_len);
        }

        let output_len = self.dims().to_len()
            * if self.output_layout() == Layout::Real {
                1
            } else {
                2
            };
        if output_len != result.len() {
            // return Err(format!("FFT plan requires that output buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", output_len).into());
            panic!("FFT plan requires that output buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", output_len);
        }

        enqueue::<T>(
            self.handle,
            direction,
            self.queue,
            buffer,
            Some(result),
            &self.wait_list,
            &mut self.dest_list,
        )
    }
}

fn enqueue<T: ClFftPrm>(
    plan: ffi::clfftPlanHandle,
    direction: Direction,
    queue: &ocl::Queue,
    buffer: &ocl::Buffer<T>,
    result: Option<&mut ocl::Buffer<T>>,
    wait_list: &Option<&ocl::EventList>,
    dest_list: &mut Option<&mut ocl::EventList>,
) -> Result<()> {
    let mut queue = queue.as_ptr();
    let mut buffer = buffer.as_core().as_ptr();
    let (wait_list_cnt, wait_list_pnt) = match *wait_list {
        None => (0, 0 as *const ocl::ffi::cl_event),
        Some(ref l) if l.len() == 0 => (0, 0 as *const ocl::ffi::cl_event),
        Some(ref l) => (l.len(), l.as_ptr() as *const ocl::ffi::cl_event),
    };
    let dest_list_pnt = match *dest_list {
        None => 0 as *mut ocl::ffi::cl_event,
        Some(ref mut l) if l.len() == 0 => 0 as *mut ocl::ffi::cl_event,
        Some(ref mut l) => unsafe { l.first_mut().unwrap().as_ptr_mut() },
    };
    let mut result = match result {
        None => 0 as ocl::ffi::cl_mem,
        Some(res) => res.as_core().as_ptr(),
    };
    let direction = translate_direction(direction);
    clfft_try!(unsafe {
        ffi::clfftEnqueueTransform(
            plan,
            direction,
            1,
            &mut queue,
            wait_list_cnt as u32,
            wait_list_pnt,
            dest_list_pnt,
            &mut buffer,
            &mut result,
            0 as ocl::ffi::cl_mem,
        )
    });
    Ok(())
}

/// Gets the native `clFFT` plan handle from a type.
pub trait AsClFftPlanHandle {
    unsafe fn as_ptr(&self) -> ffi::clfftPlanHandle;
}

impl<'a, T: ClFftPrm> AsClFftPlanHandle for FftOutOfPlacePlan<'a, T> {
    /// Returns the native clFFT plan handle.
    unsafe fn as_ptr(&self) -> ffi::clfftPlanHandle {
        self.handle
    }
}

impl<'a, T: ClFftPrm> AsClFftPlanHandle for FftInplacePlan<'a, T> {
    /// Returns the native clFFT plan handle.
    unsafe fn as_ptr(&self) -> ffi::clfftPlanHandle {
        self.handle
    }
}

/// Getters for a FFT plan.
pub trait FftPlan {
    /// Gets expected precision of each FFT.
    fn precision(&self) -> Precision;
    /// Gets the FFT dimensions.
    fn dims(&self) -> ocl::SpatialDims;
    /// the expected layouts of the input buffers.
    fn input_layout(&self) -> Layout;
    /// the expected layouts of the output buffers.
    fn output_layout(&self) -> Layout;
    /// Gets the scaling factors for FFTs.
    fn forward_scale(&self) -> f32;
    /// Gets the scaling factors for IFFTs.
    fn backward_scale(&self) -> f32;
    /// Gets the patch size.
    fn batch_size(&self) -> usize;
    /// Gets whether the input buffers are overwritten with results.
    fn result_location(&self) -> Location;
}

impl<T: AsClFftPlanHandle> FftPlan for T {
    fn precision(&self) -> Precision {
        let handle = unsafe { self.as_ptr() };
        let mut precision = ffi::clfftPrecision::CLFFT_SINGLE;
        clfft_panic!(unsafe { ffi::clfftGetPlanPrecision(handle, &mut precision) });
        translate_precision_back(precision)
    }

    fn result_location(&self) -> Location {
        let handle = unsafe { self.as_ptr() };
        let mut location = ffi::clfftResultLocation::CLFFT_INPLACE;
        clfft_panic!(unsafe { ffi::clfftGetResultLocation(handle, &mut location) });
        translate_location_back(location)
    }

    fn dims(&self) -> ocl::SpatialDims {
        let handle = unsafe { self.as_ptr() };
        let mut dim = ffi::clfftDim::CLFFT_1D;
        let mut num_dims = 0;
        clfft_panic!(unsafe { ffi::clfftGetPlanDim(handle, &mut dim, &mut num_dims) });
        match num_dims {
            1 => {
                let mut dims = [0; 1];
                clfft_panic!(unsafe { ffi::clfftGetPlanLength(handle, dim, dims.as_mut_ptr()) });
                ocl::SpatialDims::from(dims)
            }
            2 => {
                let mut dims = [0; 2];
                clfft_panic!(unsafe { ffi::clfftGetPlanLength(handle, dim, dims.as_mut_ptr()) });
                ocl::SpatialDims::from(dims)
            }
            3 => {
                let mut dims = [0; 3];
                clfft_panic!(unsafe { ffi::clfftGetPlanLength(handle, dim, dims.as_mut_ptr()) });
                ocl::SpatialDims::from(dims)
            }
            n => panic!("Unexpeced number of dimensions {}", n),
        }
    }

    fn input_layout(&self) -> Layout {
        let handle = unsafe { self.as_ptr() };
        let mut input_layout = ffi::clfftLayout::CLFFT_COMPLEX_INTERLEAVED;
        let mut output_layout = ffi::clfftLayout::CLFFT_COMPLEX_INTERLEAVED;
        clfft_panic!(unsafe { ffi::clfftGetLayout(handle, &mut input_layout, &mut output_layout) });
        translate_layout_back(input_layout)
    }

    fn output_layout(&self) -> Layout {
        let handle = unsafe { self.as_ptr() };
        let mut input_layout = ffi::clfftLayout::CLFFT_COMPLEX_INTERLEAVED;
        let mut output_layout = ffi::clfftLayout::CLFFT_COMPLEX_INTERLEAVED;
        clfft_panic!(unsafe { ffi::clfftGetLayout(handle, &mut input_layout, &mut output_layout) });
        translate_layout_back(output_layout)
    }

    fn forward_scale(&self) -> f32 {
        let handle = unsafe { self.as_ptr() };
        let mut scale = 1.0;
        clfft_panic!(unsafe {
            ffi::clfftGetPlanScale(handle, ffi::clfftDirection::CLFFT_FORWARD, &mut scale)
        });
        scale
    }

    fn backward_scale(&self) -> f32 {
        let handle = unsafe { self.as_ptr() };
        let mut scale = 1.0;
        clfft_panic!(unsafe {
            ffi::clfftGetPlanScale(handle, ffi::clfftDirection::CLFFT_BACKWARD, &mut scale)
        });
        scale
    }

    fn batch_size(&self) -> usize {
        let handle = unsafe { self.as_ptr() };
        let mut batch_size = 0;
        clfft_panic!(unsafe { ffi::clfftGetPlanBatchSize(handle, &mut batch_size) });
        batch_size
    }
}

impl<'a, T: ClFftPrm> Drop for FftInplacePlan<'a, T> {
    fn drop(&mut self) {
        if self.handle != 0 {
            let _ = unsafe { ffi::clfftDestroyPlan(&mut self.handle) };
            self.handle = 0;
        }
    }
}

impl<'a, T: ClFftPrm> Drop for FftOutOfPlacePlan<'a, T> {
    fn drop(&mut self) {
        if self.handle != 0 {
            let _ = unsafe { ffi::clfftDestroyPlan(&mut self.handle) };
            self.handle = 0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ocl::{Buffer, MemFlags};

    #[test]
    fn clfft_trivial() -> Result<()> {
        let setup = SetupData::new()?;

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
        let mut in_buffer = Buffer::builder()
            .queue(queue.clone())
            .flags(MemFlags::new().read_write().copy_host_ptr())
            .len(source.len())
            .copy_host_slice(&source)
            .build()?;

        let mut res_buffer = Buffer::<f64>::builder()
            .queue(queue.clone())
            .flags(MemFlags::new().write_only())
            .len(2 * source.len())
            .build()
            .expect("Failed to create GPU result buffer");

        // Make a plan
        let mut plan = builder::<f64>(&setup)
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

        let setup = SetupData::new()?;

        let ctx = ocl::Context::builder().build()?;
        let device = ocl::Device::first(ocl::Platform::first()?)?;
        let queue = ocl::Queue::new(&ctx, device, None)?;

        let x = vec![0f32; 2 * N];
        let mut bufx = ocl::Buffer::builder()
            .queue(queue.clone())
            .flags(ocl::MemFlags::new().read_write())
            .len(2 * N)
            .copy_host_slice(&x)
            .build()?;

        bufx.cmd().write(&x).enq()?;

        // // Make a plan
        let mut plan = builder::<f32>(&setup)
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
}
