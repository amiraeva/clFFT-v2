use std::ptr;

// #[link(name = "clFFT", kind = "dylib")]
pub mod ffi;

const DIM: ffi::clfftDim = ffi::clfftDim::CLFFT_1D;

fn init_setup_data() -> ffi::clfftSetupData {
    let mut sdata = ffi::clfftSetupData::default();

    sdata.major = ffi::clfftVersionMajor;
    sdata.minor = ffi::clfftVersionMinor;
    sdata.patch = ffi::clfftVersionPatch;
    sdata.debugFlags = 0;

    sdata
}
#[derive(Debug)]
enum Error {
    ffi(ffi::clfftStatus),
    ocl(ocl::Error),
}

impl From<ffi::clfftStatus> for Result<(), Error> {
    fn from(val: ffi::clfftStatus) -> Self {
        match val {
            ffi::clfftStatus::CLFFT_SUCCESS => Ok(()),
            _ => Err(Error::ffi(val)),
        }
    }
}

impl From<ocl::Error> for Error {
    fn from(val: ocl::Error) -> Self {
        Self::ocl(val)
    }
}

struct SetupData(ffi::clfftSetupData);

impl SetupData {
    fn new() -> Result<Self, Error> {
        let sdata = init_setup_data();
        let err = unsafe { ffi::clfftSetup(&sdata) };
        Result::from(err)?;
        Ok(Self(sdata))
    }
}
struct PlanHandle(ffi::clfftPlanHandle);

impl PlanHandle {
    fn new() -> Self {
        Self(ffi::clfftPlanHandle::default())
    }

    fn create_default(&mut self, ctx: &ocl::Context, len: usize) -> Result<&mut Self, Error> {
        let err = unsafe { ffi::clfftCreateDefaultPlan(&mut self.0, ctx.as_ptr(), DIM, &len) };
        Result::from(err)?;
        Ok(self)
    }

    fn set_precision(&mut self, precision: ffi::clfftPrecision) -> Result<&mut Self, Error> {
        let err = unsafe { ffi::clfftSetPlanPrecision(self.0, precision) };
        Result::from(err)?;
        Ok(self)
    }

    fn set_layout(
        &mut self,
        iLayout: ffi::clfftLayout,
        oLayout: ffi::clfftLayout,
    ) -> Result<&mut Self, Error> {
        let err = unsafe { ffi::clfftSetLayout(self.0, iLayout, oLayout) };
        Result::from(err)?;
        Ok(self)
    }

    fn set_result_location(&mut self, place: ffi::clfftResultLocation) -> Result<&mut Self, Error> {
        let err = unsafe { ffi::clfftSetResultLocation(self.0, place) };
        Result::from(err)?;
        Ok(self)
    }

    fn bake_plan(&mut self, queue: &mut ocl::Queue) -> Result<&mut Self, Error> {
        let err =
            unsafe { ffi::clfftBakePlan(self.0, 1, &mut queue.as_ptr(), None, ptr::null_mut()) };
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clfft_trivial() -> Result<(), Error> {
        const N: usize = 16;

        let setup = SetupData::new();

        let mut phandle = ffi::clfftPlanHandle::default();
        let mut pro_que = ocl::ProQue::builder().build()?;

        let mut phandle = PlanHandle::new();

        phandle
            .create_default(pro_que.context(), N)?
            .set_precision(ffi::clfftPrecision::CLFFT_SINGLE_FAST)?
            .set_layout(
                ffi::clfftLayout::CLFFT_REAL,
                ffi::clfftLayout::CLFFT_HERMITIAN_INTERLEAVED,
            )?
            .set_result_location(ffi::clfftResultLocation::CLFFT_INPLACE);

        // err = ffi::clfftCreateDefaultPlan(&mut phandle, ctx.as_ptr(), DIM, &N);
        // err = ffi::clfftSetPlanPrecision(phandle, ffi::clfftPrecision::CLFFT_SINGLE_FAST);

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
