use std::marker::PhantomData;

pub mod ffi;
mod util;

macro_rules! clfft_panic {
    ( $result_expr: expr) => {
        let result = $result_expr;
        if result != ffi::clfftStatus::CLFFT_SUCCESS {
            panic!("Unexpeced error in foreign library: {:?}", result);
        }
    };
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

pub struct SetupData(ffi::clfftSetupData);

impl SetupData {
    pub fn new() -> Result<Self> {
        let sdata = util::init_setup_data();
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

/// Specify the expected precision of each FFT.
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Precision {
    Precise,
    Fast,
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

/// Specify the expected direction of each FFT, time or the frequency domains
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Direction {
    Forward,
    Backward,
}

/// pecify wheter the input buffers are overwritten with results
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Location {
    Inplace,
    OutOfPlace,
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
    ) -> Result<FFTPlan<'a, T, InPlace>> {
        let handle = util::bake_plan::<T>(self, queue, ctx, Location::Inplace)?;
        Ok(FFTPlan {
            handle: handle,
            queue,
            data_type: std::marker::PhantomData,
            wait_list: None,
            dest_list: None,
            _location: InPlace,
        })
    }

    /// Creates a plan for an out of place FFT.
    pub fn bake_out_of_place_plan<'a>(
        &mut self,
        queue: &'a ocl::Queue,
        ctx: &ocl::Context,
    ) -> Result<FFTPlan<'a, T, OutOfPlace>> {
        let handle = util::bake_plan::<T>(self, queue, ctx, Location::OutOfPlace)?;
        Ok(FFTPlan {
            handle: handle,
            queue,
            data_type: std::marker::PhantomData,
            wait_list: None,
            dest_list: None,
            _location: OutOfPlace,
        })
    }
}

pub struct InPlace;
pub struct OutOfPlace;

pub trait Placeness {}

impl Placeness for InPlace {}
impl Placeness for OutOfPlace {}

pub struct FFTPlan<'a, T: ClFftPrm, L: Placeness> {
    handle: ffi::clfftPlanHandle,
    queue: &'a ocl::Queue,
    data_type: std::marker::PhantomData<T>,
    wait_list: Option<&'a ocl::EventList>,
    dest_list: Option<&'a mut ocl::EventList>,
    _location: L,
}

impl<'a, T: ClFftPrm, L: Placeness> FFTPlan<'a, T, L> {
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

    unsafe fn as_ptr(&self) -> ffi::clfftPlanHandle {
        self.handle
    }

    pub fn precision(&self) -> Precision {
        let handle = unsafe { self.as_ptr() };
        let mut precision = ffi::clfftPrecision::CLFFT_SINGLE;
        clfft_panic!(unsafe { ffi::clfftGetPlanPrecision(handle, &mut precision) });
        util::translate_precision_back(precision)
    }

    pub fn result_location(&self) -> Location {
        let handle = unsafe { self.as_ptr() };
        let mut location = ffi::clfftResultLocation::CLFFT_INPLACE;
        clfft_panic!(unsafe { ffi::clfftGetResultLocation(handle, &mut location) });
        util::translate_location_back(location)
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
        util::translate_layout_back(input_layout)
    }

    fn output_layout(&self) -> Layout {
        let handle = unsafe { self.as_ptr() };
        let mut input_layout = ffi::clfftLayout::CLFFT_COMPLEX_INTERLEAVED;
        let mut output_layout = ffi::clfftLayout::CLFFT_COMPLEX_INTERLEAVED;
        clfft_panic!(unsafe { ffi::clfftGetLayout(handle, &mut input_layout, &mut output_layout) });
        util::translate_layout_back(output_layout)
    }

    pub fn forward_scale(&self) -> f32 {
        let handle = unsafe { self.as_ptr() };
        let mut scale = 1.0;
        clfft_panic!(unsafe {
            ffi::clfftGetPlanScale(handle, ffi::clfftDirection::CLFFT_FORWARD, &mut scale)
        });
        scale
    }

    pub fn backward_scale(&self) -> f32 {
        let handle = unsafe { self.as_ptr() };
        let mut scale = 1.0;
        clfft_panic!(unsafe {
            ffi::clfftGetPlanScale(handle, ffi::clfftDirection::CLFFT_BACKWARD, &mut scale)
        });
        scale
    }

    pub fn batch_size(&self) -> usize {
        let handle = unsafe { self.as_ptr() };
        let mut batch_size = 0;
        clfft_panic!(unsafe { ffi::clfftGetPlanBatchSize(handle, &mut batch_size) });
        batch_size
    }

    fn enq_inner(
        &mut self,
        direction: Direction,
        buffer: &ocl::Buffer<T>,
        result: Option<&mut ocl::Buffer<T>>,
    ) -> Result<()> {
        util::enqueue(
            self.handle,
            direction,
            self.queue,
            buffer,
            result,
            &self.wait_list,
            &mut self.dest_list,
        )
    }
}

impl<'a, T: ClFftPrm> FFTPlan<'a, T, InPlace> {
    /// Enqueues the FFT so that it gets performed on the device.
    pub fn enq(&mut self, direction: Direction, buffer: &mut ocl::Buffer<T>) -> Result<()> {
        util::validate_buffer_len(buffer, self.input_layout(), self.dims());
        self.enq_inner(direction, buffer, None)
    }
}

impl<'a, T: ClFftPrm> FFTPlan<'a, T, OutOfPlace> {
    /// Enqueues the FFT so that it gets performed on the device.
    pub fn enq(
        &mut self,
        direction: Direction,
        buffer: &ocl::Buffer<T>,
        result: &mut ocl::Buffer<T>,
    ) -> Result<()> {
        util::validate_buffer_len(buffer, self.input_layout(), self.dims());
        util::validate_buffer_len(result, self.output_layout(), self.dims());

        self.enq_inner(direction, buffer, Some(result))
    }
}

impl<'a, T: ClFftPrm, L: Placeness> Drop for FFTPlan<'a, T, L> {
    fn drop(&mut self) {
        if self.handle != 0 {
            let _ = unsafe { ffi::clfftDestroyPlan(&mut self.handle) };
            self.handle = 0;
        }
    }
}

#[cfg(test)]
mod tests;
