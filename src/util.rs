use super::*;

use ffi::*;
use ocl::ffi::{cl_event, cl_mem};

macro_rules! clfft_try {
    ( $result_expr: expr) => {
        let err = $result_expr;
        Result::from(err)?;
    };
}

pub(crate) fn init_setup_data() -> clfftSetupData {
    let mut sdata = clfftSetupData::default();

    sdata.major = clfftVersionMajor;
    sdata.minor = clfftVersionMinor;
    sdata.patch = clfftVersionPatch;
    sdata.debugFlags = 0;

    sdata
}

fn translate_to_fft_dim(dims: ocl::SpatialDims) -> clfftDim {
    match dims.dim_count() {
        1 => clfftDim::CLFFT_1D,
        2 => clfftDim::CLFFT_2D,
        3 => clfftDim::CLFFT_3D,
        n => panic!("Number of dimensions must be 1, 2 or 3, but it is {}", n),
    }
}

fn translate_precision<T: ClFftPrm>(precision: Precision) -> clfftPrecision {
    let is_f64 = T::is_dbl_precision();
    match precision {
        Precision::Precise => {
            if is_f64 {
                clfftPrecision::CLFFT_DOUBLE
            } else {
                clfftPrecision::CLFFT_SINGLE
            }
        }
        Precision::Fast => {
            if is_f64 {
                clfftPrecision::CLFFT_DOUBLE_FAST
            } else {
                clfftPrecision::CLFFT_SINGLE_FAST
            }
        }
    }
}

pub(crate) fn translate_precision_back(precision: clfftPrecision) -> Precision {
    match precision {
        clfftPrecision::CLFFT_SINGLE | clfftPrecision::CLFFT_DOUBLE => Precision::Precise,
        clfftPrecision::CLFFT_SINGLE_FAST | clfftPrecision::CLFFT_DOUBLE_FAST => Precision::Fast,
        clfftPrecision::ENDPRECISION => unreachable!("ENDPRECISION should never be returned"),
    }
}

fn translate_layout(layout: Layout) -> clfftLayout {
    match layout {
        Layout::ComplexInterleaved => clfftLayout::CLFFT_COMPLEX_INTERLEAVED,
        Layout::ComplexPlanar => clfftLayout::CLFFT_COMPLEX_PLANAR,
        Layout::HermitianInterleaved => clfftLayout::CLFFT_HERMITIAN_INTERLEAVED,
        Layout::HermitianPlanar => clfftLayout::CLFFT_HERMITIAN_PLANAR,
        Layout::Real => clfftLayout::CLFFT_REAL,
    }
}

pub(crate) fn translate_layout_back(layout: clfftLayout) -> Layout {
    match layout {
        clfftLayout::CLFFT_COMPLEX_INTERLEAVED => Layout::ComplexInterleaved,
        clfftLayout::CLFFT_COMPLEX_PLANAR => Layout::ComplexPlanar,
        clfftLayout::CLFFT_HERMITIAN_INTERLEAVED => Layout::HermitianInterleaved,
        clfftLayout::CLFFT_HERMITIAN_PLANAR => Layout::HermitianPlanar,
        clfftLayout::CLFFT_REAL => Layout::Real,
        clfftLayout::ENDLAYOUT => panic!("ENDLAYOUT should never be returned"),
    }
}

fn translate_direction(direction: Direction) -> clfftDirection {
    match direction {
        Direction::Forward => clfftDirection::CLFFT_FORWARD,
        Direction::Backward => clfftDirection::CLFFT_BACKWARD,
    }
}

fn translate_location(location: Location) -> clfftResultLocation {
    match location {
        Location::Inplace => clfftResultLocation::CLFFT_INPLACE,
        Location::OutOfPlace => clfftResultLocation::CLFFT_OUTOFPLACE,
    }
}

pub(crate) fn translate_location_back(location: clfftResultLocation) -> Location {
    match location {
        clfftResultLocation::CLFFT_INPLACE => Location::Inplace,
        clfftResultLocation::CLFFT_OUTOFPLACE => Location::OutOfPlace,
        clfftResultLocation::ENDPLACE => panic!("ENDPLACE should never be returned"),
    }
}

pub(crate) fn bake_plan<T: ClFftPrm>(
    builder: &mut FftPlanBuilder<T>,
    queue: &ocl::Queue,
    ctx: &ocl::Context,
    location: Location,
) -> Result<clfftPlanHandle> {
    let context = ctx.as_ptr();
    let dims = builder.dims.unwrap();

    let dim = util::translate_to_fft_dim(dims);
    let lengths = dims
        .to_lens()
        .map_err(|_| Error::UnspecifiedDimensionError)?;

    let mut plan: clfftPlanHandle = 0;
    clfft_try!(unsafe { clfftCreateDefaultPlan(&mut plan, context, dim, lengths.as_ptr()) });
    let precision = translate_precision::<T>(builder.precision);
    clfft_try!(unsafe { clfftSetPlanPrecision(plan, precision) });
    let input_layout = translate_layout(builder.input_layout);
    let output_layout = translate_layout(builder.output_layout);
    clfft_try!(unsafe { clfftSetLayout(plan, input_layout, output_layout) });
    match builder.forward_scale {
        None => (),
        Some(s) => {
            clfft_try!(unsafe { clfftSetPlanScale(plan, clfftDirection::CLFFT_FORWARD, s) });
        }
    }
    match builder.backward_scale {
        None => (),
        Some(s) => {
            clfft_try!(unsafe { clfftSetPlanScale(plan, clfftDirection::CLFFT_BACKWARD, s) });
        }
    }
    let location = translate_location(location);
    clfft_try!(unsafe { clfftSetResultLocation(plan, location) });
    match builder.batch_size {
        None => (), // Use default
        Some(s) => {
            clfft_try!(unsafe { clfftSetPlanBatchSize(plan, s) });
        }
    }

    let mut queue = queue.as_ptr();
    clfft_try!(unsafe {
        clfftBakePlan(plan, 1, &mut queue, None, 0 as *mut ::std::os::raw::c_void)
    });
    Ok(plan)
}

pub(crate) fn enqueue<T: ClFftPrm>(
    plan: clfftPlanHandle,
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
        None => (0, 0 as *const cl_event),
        Some(ref l) if l.len() == 0 => (0, 0 as *const cl_event),
        Some(ref l) => (l.len(), l.as_ptr() as *const cl_event),
    };
    let dest_list_pnt = match *dest_list {
        None => 0 as *mut cl_event,
        Some(ref mut l) if l.len() == 0 => 0 as *mut cl_event,
        Some(ref mut l) => unsafe { l.first_mut().unwrap().as_ptr_mut() },
    };
    let mut result = match result {
        None => 0 as cl_mem,
        Some(res) => res.as_core().as_ptr(),
    };
    let direction = translate_direction(direction);
    clfft_try!(unsafe {
        clfftEnqueueTransform(
            plan,
            direction,
            1,
            &mut queue,
            wait_list_cnt as u32,
            wait_list_pnt,
            dest_list_pnt,
            &mut buffer,
            &mut result,
            0 as cl_mem,
        )
    });
    Ok(())
}

pub(crate) fn validate_buffer_len<T: ocl::OclPrm>(buffer: &ocl::Buffer<T>, buf_layout: Layout, dims: ocl::SpatialDims) {
	let factor = if buf_layout == Layout::Real {
		1
	} else {
		2
	};

	let expected_len = dims.to_len() * factor;

	// TODO get more meaningful error messages
	// if input_len != buffer.len() {
	//     // return Err(format!("FFT plan requires that input buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", input_len).into());
	//     // TODO: convert to error
	//     panic!("FFT plan requires that input buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", input_len);
	// }

	// if output_len != result.len() {
	//     // return Err(format!("FFT plan requires that output buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", output_len).into());
	//     // TODO: convert to error
	//     panic!("FFT plan requires that output buffer must have a size of {}. Is there a dimension mismatch between real and complex numbers?", output_len);
	// }

	if expected_len != buffer.len() {
		// TODO: convert to error
		panic!("Input or output buffer length mismatch: expected {}, found: {}", expected_len, buffer.len());
	}
}
