#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use cl_sys::{cl_context, cl_mem, cl_ulong, cl_float, cl_command_queue, cl_uint, cl_event};

pub const clfftVersionMajor: cl_uint = 2;
pub const clfftVersionMinor: cl_uint = 12;
pub const clfftVersionPatch: cl_uint = 2;

#[repr(i32)]
#[doc = "  @brief clfft error codes definition(incorporating OpenCL error definitions)"]
#[doc = ""]
#[doc = "   This enumeration is a superset of the OpenCL error codes.  For example, CL_OUT_OF_HOST_MEMORY,"]
#[doc = "   which is defined in cl.h is aliased as CLFFT_OUT_OF_HOST_MEMORY.  The set of basic OpenCL"]
#[doc = "   error codes is extended to add extra values specific to the clfft package."]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum clfftStatus {
    CLFFT_INVALID_GLOBAL_WORK_SIZE = -63,
    CLFFT_INVALID_MIP_LEVEL = -62,
    CLFFT_INVALID_BUFFER_SIZE = -61,
    CLFFT_INVALID_GL_OBJECT = -60,
    CLFFT_INVALID_OPERATION = -59,
    CLFFT_INVALID_EVENT = -58,
    CLFFT_INVALID_EVENT_WAIT_LIST = -57,
    CLFFT_INVALID_GLOBAL_OFFSET = -56,
    CLFFT_INVALID_WORK_ITEM_SIZE = -55,
    CLFFT_INVALID_WORK_GROUP_SIZE = -54,
    CLFFT_INVALID_WORK_DIMENSION = -53,
    CLFFT_INVALID_KERNEL_ARGS = -52,
    CLFFT_INVALID_ARG_SIZE = -51,
    CLFFT_INVALID_ARG_VALUE = -50,
    CLFFT_INVALID_ARG_INDEX = -49,
    CLFFT_INVALID_KERNEL = -48,
    CLFFT_INVALID_KERNEL_DEFINITION = -47,
    CLFFT_INVALID_KERNEL_NAME = -46,
    CLFFT_INVALID_PROGRAM_EXECUTABLE = -45,
    CLFFT_INVALID_PROGRAM = -44,
    CLFFT_INVALID_BUILD_OPTIONS = -43,
    CLFFT_INVALID_BINARY = -42,
    CLFFT_INVALID_SAMPLER = -41,
    CLFFT_INVALID_IMAGE_SIZE = -40,
    CLFFT_INVALID_IMAGE_FORMAT_DESCRIPTOR = -39,
    CLFFT_INVALID_MEM_OBJECT = -38,
    CLFFT_INVALID_HOST_PTR = -37,
    CLFFT_INVALID_COMMAND_QUEUE = -36,
    CLFFT_INVALID_QUEUE_PROPERTIES = -35,
    CLFFT_INVALID_CONTEXT = -34,
    CLFFT_INVALID_DEVICE = -33,
    CLFFT_INVALID_PLATFORM = -32,
    CLFFT_INVALID_DEVICE_TYPE = -31,
    CLFFT_INVALID_VALUE = -30,
    CLFFT_MAP_FAILURE = -12,
    CLFFT_BUILD_PROGRAM_FAILURE = -11,
    CLFFT_IMAGE_FORMAT_NOT_SUPPORTED = -10,
    CLFFT_IMAGE_FORMAT_MISMATCH = -9,
    CLFFT_MEM_COPY_OVERLAP = -8,
    CLFFT_PROFILING_INFO_NOT_AVAILABLE = -7,
    CLFFT_OUT_OF_HOST_MEMORY = -6,
    CLFFT_OUT_OF_RESOURCES = -5,
    CLFFT_MEM_OBJECT_ALLOCATION_FAILURE = -4,
    CLFFT_COMPILER_NOT_AVAILABLE = -3,
    CLFFT_DEVICE_NOT_AVAILABLE = -2,
    CLFFT_DEVICE_NOT_FOUND = -1,
    CLFFT_SUCCESS = 0,
    #[doc = "< Bugcheck."]
    CLFFT_BUGCHECK = 4096,
    #[doc = "< Functionality is not implemented yet."]
    CLFFT_NOTIMPLEMENTED = 4097,
    #[doc = "< Transposed functionality is not implemented for this transformation."]
    CLFFT_TRANSPOSED_NOTIMPLEMENTED = 4098,
    #[doc = "< Tried to open an existing file on the host system, but failed."]
    CLFFT_FILE_NOT_FOUND = 4099,
    #[doc = "< Tried to create a file on the host system, but failed."]
    CLFFT_FILE_CREATE_FAILURE = 4100,
    #[doc = "< Version conflict between client and library."]
    CLFFT_VERSION_MISMATCH = 4101,
    #[doc = "< Requested plan could not be found."]
    CLFFT_INVALID_PLAN = 4102,
    #[doc = "< Double precision not supported on this device."]
    CLFFT_DEVICE_NO_DOUBLE = 4103,
    #[doc = "< Attempt to run on a device using a plan baked for a different device."]
    CLFFT_DEVICE_MISMATCH = 4104,
    CLFFT_ENDSTATUS = 4105,
}

#[repr(u32)]
#[doc = "  @brief The dimension of the input and output buffers that is fed into all FFT transforms"]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum clfftDim {
    #[doc = "< 1 Dimensional FFT transform (default)."]
    CLFFT_1D = 1,
    #[doc = "< 2 Dimensional FFT transform."]
    CLFFT_2D = 2,
    #[doc = "< 3 Dimensional FFT transform."]
    CLFFT_3D = 3,
    #[doc = "< The last value of the enum, and marks the length of clfftDim."]
    ENDDIMENSION = 4,
}
#[repr(u32)]
#[doc = "  @brief Specify the expected layouts of the buffers"]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum clfftLayout_ {
    #[doc = "< An array of complex numbers, with real and imaginary components together (default)."]
    CLFFT_COMPLEX_INTERLEAVED = 1,
    #[doc = "< Separate arrays of real components and imaginary components."]
    CLFFT_COMPLEX_PLANAR = 2,
    #[doc = "< Compressed form of complex numbers; complex-conjugates are not stored, real and imaginary components are stored in the same array."]
    CLFFT_HERMITIAN_INTERLEAVED = 3,
    #[doc = "< Compressed form of complex numbers; complex-conjugates are not stored, real and imaginary components are stored in separate arrays."]
    CLFFT_HERMITIAN_PLANAR = 4,
    #[doc = "< An array of real numbers, with no corresponding imaginary components."]
    CLFFT_REAL = 5,
    #[doc = "< The last value of the enum, and marks the length of clfftLayout."]
    ENDLAYOUT = 6,
}
#[doc = "  @brief Specify the expected layouts of the buffers"]
pub use self::clfftLayout_ as clfftLayout;
#[repr(u32)]
#[doc = "  @brief Specify the expected precision of each FFT."]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum clfftPrecision {
    #[doc = "< An array of complex numbers, with real and imaginary components saved as floats (default)."]
    CLFFT_SINGLE = 1,
    #[doc = "< An array of complex numbers, with real and imaginary components saved as doubles."]
    CLFFT_DOUBLE = 2,
    #[doc = "< Faster implementation preferred."]
    CLFFT_SINGLE_FAST = 3,
    #[doc = "< Faster implementation preferred."]
    CLFFT_DOUBLE_FAST = 4,
    #[doc = "< The last value of the enum, and marks the length of clfftPrecision."]
    ENDPRECISION = 5,
}

impl clfftDirection_ {
    pub const CLFFT_MINUS: clfftDirection_ = clfftDirection_::CLFFT_FORWARD;
}
impl clfftDirection_ {
    pub const CLFFT_PLUS: clfftDirection_ = clfftDirection_::CLFFT_BACKWARD;
}
#[repr(i32)]
#[doc = "  @brief Specify the expected direction of each FFT, time or the frequency domains"]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum clfftDirection_ {
    #[doc = "< FFT transform from time to frequency domain."]
    CLFFT_FORWARD = -1,
    #[doc = "< FFT transform from frequency to time domain."]
    CLFFT_BACKWARD = 1,
    #[doc = "< The last value of the enum, and marks the length of clfftDirection."]
    ENDDIRECTION = 2,
}
#[doc = "  @brief Specify the expected direction of each FFT, time or the frequency domains"]
pub use self::clfftDirection_ as clfftDirection;
#[repr(u32)]
#[doc = "  @brief Specify wheter the input buffers are overwritten with results"]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum clfftResultLocation_ {
    #[doc = "< Input and output buffers are the same (default)."]
    CLFFT_INPLACE = 1,
    #[doc = "< Input and output buffers are separate."]
    CLFFT_OUTOFPLACE = 2,
    #[doc = "< The last value of the enum, and marks the length of clfftPlaceness."]
    ENDPLACE = 3,
}
#[doc = "  @brief Specify wheter the input buffers are overwritten with results"]
pub use self::clfftResultLocation_ as clfftResultLocation;
#[repr(u32)]
#[doc = " @brief Determines whether the result is returned in original order. It is valid only for"]
#[doc = "dimensions greater than 1."]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum clfftResultTransposed_ {
    #[doc = "< The result is returned in the original order (default)"]
    CLFFT_NOTRANSPOSE = 1,
    #[doc = "< The result is transposed where transpose kernel is supported (possibly faster)"]
    CLFFT_TRANSPOSED = 2,
    #[doc = "< The last value of the enum, and marks the length of clfftResultTransposed"]
    ENDTRANSPOSED = 3,
}
#[doc = " @brief Determines whether the result is returned in original order. It is valid only for"]
#[doc = "dimensions greater than 1."]
pub use self::clfftResultTransposed_ as clfftResultTransposed;
#[doc = " @brief Data structure that can be passed to clfftSetup() to control the behavior of the FFT runtime"]
#[doc = "  @details This structure contains values that can be initialized before instantiation of the FFT runtime"]
#[doc = "  with ::clfftSetup().  To initialize this structure, pass a pointer to a user struct to ::clfftInitSetupData( ),"]
#[doc = "  which clears the structure and sets the version member variables to the current values."]
#[repr(C)]
#[derive(Debug, Default, Copy, Clone)]
pub struct clfftSetupData {
    #[doc = "< Major version number of the project; signifies possible major API changes."]
    pub major: cl_uint,
    #[doc = "< Minor version number of the project; minor API changes that can break backward compatibility."]
    pub minor: cl_uint,
    #[doc = "< Patch version number of the project; always incrementing number, signifies change over time."]
    pub patch: cl_uint,
    #[doc = " \tBitwise flags that control the behavior of library debug logic."]
    pub debugFlags: cl_ulong,
}
#[test]
fn bindgen_test_layout_clfftSetupData() {
    assert_eq!(
        ::std::mem::size_of::<clfftSetupData>(),
        24usize,
        concat!("Size of: ", stringify!(clfftSetupData))
    );
    assert_eq!(
        ::std::mem::align_of::<clfftSetupData>(),
        8usize,
        concat!("Alignment of ", stringify!(clfftSetupData))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<clfftSetupData>())).major as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(clfftSetupData),
            "::",
            stringify!(major)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<clfftSetupData>())).minor as *const _ as usize },
        4usize,
        concat!(
            "Offset of field: ",
            stringify!(clfftSetupData),
            "::",
            stringify!(minor)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<clfftSetupData>())).patch as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(clfftSetupData),
            "::",
            stringify!(patch)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<clfftSetupData>())).debugFlags as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(clfftSetupData),
            "::",
            stringify!(debugFlags)
        )
    );
}

#[repr(u32)]
#[doc = " @brief Type of Callback function."]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum clfftCallbackType_ {
    #[doc = "< Callback function is invoked only once for every point of input at the beginning of FFT transform."]
    PRECALLBACK = 0,
    #[doc = "< Callback function is invoked only once for every point of output at the end of FFT transform."]
    POSTCALLBACK = 1,
}
#[doc = " @brief Type of Callback function."]
pub use self::clfftCallbackType_ as clfftCallbackType;
#[doc = "  @brief An abstract handle to the object that represents the state of the FFT(s)"]
pub type clfftPlanHandle = usize;
extern "C" {
    #[doc = " @brief Initialize the internal FFT resources."]
    #[doc = "  @details The internal resources include FFT implementation caches kernels, programs, and buffers."]
    #[doc = "  @param[in] setupData Data structure that is passed into the setup routine to control FFT generation behavior"]
    #[doc = " \tand debug functionality"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftSetup(setupData: *const clfftSetupData) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Release all internal resources."]
    #[doc = "  @details Called when client is done with the FFT library, allowing the library to destroy all resources it has cached"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftTeardown() -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Query the FFT library for version information"]
    #[doc = "  @details Returns the major, minor and patch version numbers associated with the FFT library"]
    #[doc = "  @param[out] major Major functionality change"]
    #[doc = "  @param[out] minor Minor functionality change"]
    #[doc = "  @param[out] patch Bug fixes, documentation changes, no new features introduced"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftGetVersion(
        major: *mut cl_uint,
        minor: *mut cl_uint,
        patch: *mut cl_uint,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Create a plan object initialized entirely with default values."]
    #[doc = "  @details A plan is a repository of state for calculating FFT's.  Allows the runtime to pre-calculate kernels, programs"]
    #[doc = " \tand buffers and associate them with buffers of specified dimensions."]
    #[doc = "  @param[out] plHandle Handle to the newly created plan"]
    #[doc = "  @param[in] context Client is responsible for providing an OpenCL context for the plan"]
    #[doc = "  @param[in] dim Dimensionality of the FFT transform; describes how many elements are in the array"]
    #[doc = "  @param[in] clLengths An array of length of size 'dim';  each array value describes the length of each dimension"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftCreateDefaultPlan(
        plHandle: *mut clfftPlanHandle,
        context: cl_context,
        dim: clfftDim,
        clLengths: *const usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Create a copy of an existing plan."]
    #[doc = "  @details This API allows a client to create a new plan based upon an existing plan.  This function can be used to"]
    #[doc = "  quickly create plans that are similar, but may differ slightly."]
    #[doc = "  @param[out] out_plHandle Handle to the newly created plan that is based on in_plHandle"]
    #[doc = "  @param[in] new_context Client is responsible for providing a new context for the new plan"]
    #[doc = "  @param[in] in_plHandle Handle to a previously created plan that is to be copied"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftCopyPlan(
        out_plHandle: *mut clfftPlanHandle,
        new_context: cl_context,
        in_plHandle: clfftPlanHandle,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Prepare the plan for execution."]
    #[doc = "  @details After all plan parameters are set, the client has the option of 'baking' the plan, which informs the runtime that"]
    #[doc = "  no more change to the parameters of the plan is expected, and the OpenCL kernels can be compiled.  This optional function"]
    #[doc = "  allows the client application to perform the OpenCL kernel compilation when the application is initialized instead of during the first"]
    #[doc = "  execution."]
    #[doc = "  At this point, the clfft runtime applies all implimented optimizations, including"]
    #[doc = "  running kernel experiments on the devices in the plan context."]
    #[doc = "  <p>  This function takes a long time to execute. If a plan is not baked before being executed,"]
    #[doc = "  the first call to clfftEnqueueTransform takes a long time to execute."]
    #[doc = "  <p>  If any significant parameter of a plan is changed after the plan is baked (by a subsequent call to any one of"]
    #[doc = "  the functions that has the prefix \"clfftSetPlan\"), it is not considered an error.  Instead, the plan reverts back to"]
    #[doc = "  the unbaked state, discarding the benefits of the baking operation."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] numQueues Number of command queues in commQueueFFT; 0 is a valid value, in which case the client does not want"]
    #[doc = " \tthe runtime to run load experiments and only pre-calculate state information"]
    #[doc = "  @param[in] commQueueFFT An array of cl_command_queues created by the client; the command queues must be a proper subset of"]
    #[doc = " \tthe devices included in the plan context"]
    #[doc = "  @param[in] pfn_notify A function pointer to a notification routine. The notification routine is a callback function that"]
    #[doc = "  an application can register and is called when the program executable is built (successfully or unsuccessfully)."]
    #[doc = "  Currently, this parameter MUST be NULL or nullptr."]
    #[doc = "  @param[in] user_data Passed as an argument when pfn_notify is called."]
    #[doc = "  Currently, this parameter MUST be NULL or nullptr."]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftBakePlan(
        plHandle: clfftPlanHandle,
        numQueues: cl_uint,
        commQueueFFT: *mut cl_command_queue,
        pfn_notify: ::std::option::Option<
            unsafe extern "C" fn(plHandle: clfftPlanHandle, user_data: *mut ::std::os::raw::c_void),
        >,
        user_data: *mut ::std::os::raw::c_void,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Release the resources of a plan."]
    #[doc = "  @details A plan may include resources, such as kernels, programs, and buffers that consume memory.  When a plan"]
    #[doc = "  is no more needed, the client must release the plan."]
    #[doc = "  @param[in,out] plHandle Handle to a previously created plan"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftDestroyPlan(plHandle: *mut clfftPlanHandle) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the OpenCL context of a previously created plan."]
    #[doc = "  @details The user must pass a reference to a cl_context variable, which is modified to point to a"]
    #[doc = "  context set in the specified plan."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] context Reference to the user allocated cl_context, which points to context set in the plan"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftGetPlanContext(plHandle: clfftPlanHandle, context: *mut cl_context) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the floating point precision of the FFT data"]
    #[doc = "  @details The user must pass a reference to a clfftPrecision variable, which is set to the"]
    #[doc = "  precision of the FFT complex data in the plan."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] precision Reference to the user clfftPrecision enum"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftGetPlanPrecision(
        plHandle: clfftPlanHandle,
        precision: *mut clfftPrecision,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the floating point precision of the FFT data"]
    #[doc = "  @details Sets the floating point precision of the FFT complex data in the plan."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] precision Reference to the user clfftPrecision enum"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftSetPlanPrecision(
        plHandle: clfftPlanHandle,
        precision: clfftPrecision,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the scaling factor that is applied to the FFT data"]
    #[doc = "  @details The user must pass a reference to a cl_float variable, which is set to the"]
    #[doc = "  floating point scaling factor that is multiplied across the FFT data."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dir Direction of the applied scaling factor"]
    #[doc = "  @param[out] scale Reference to the user cl_float variable"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftGetPlanScale(
        plHandle: clfftPlanHandle,
        dir: clfftDirection,
        scale: *mut cl_float,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the scaling factor that is applied to the FFT data"]
    #[doc = "  @details Sets the floating point scaling factor that is"]
    #[doc = "  multiplied across the FFT data."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dir Direction of the applied scaling factor"]
    #[doc = "  @param[in] scale Reference to the user cl_float variable"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftSetPlanScale(
        plHandle: clfftPlanHandle,
        dir: clfftDirection,
        scale: cl_float,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the number of discrete arrays that the plan can concurrently handle"]
    #[doc = "  @details The user must pass a reference to a cl_uint variable, which is set to the"]
    #[doc = "  number of discrete arrays (1D or 2D) that is batched together for the plan"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] batchSize Number of discrete FFTs performed"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftGetPlanBatchSize(plHandle: clfftPlanHandle, batchSize: *mut usize) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the number of discrete arrays that the plan can concurrently handle"]
    #[doc = "  @details Sets the plan property which sets the number of discrete arrays (1D or 2D)"]
    #[doc = "  that is batched together for the plan"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] batchSize Number of discrete FFTs performed"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftSetPlanBatchSize(plHandle: clfftPlanHandle, batchSize: usize) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the dimensionality of the data that is transformed"]
    #[doc = "  @details Queries a plan object and retrieves the value of the dimensionality that the plan is set for.  A size is returned to"]
    #[doc = "  help the client allocate sufficient storage to hold the dimensions in a further call to clfftGetPlanLength"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] dim The dimensionality of the FFT to be transformed"]
    #[doc = "  @param[out] size Value to allocate an array to hold the FFT dimensions."]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftGetPlanDim(
        plHandle: clfftPlanHandle,
        dim: *mut clfftDim,
        size: *mut cl_uint,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the dimensionality of the data that is transformed"]
    #[doc = "  @details Set the dimensionality of the data that is transformed by the plan"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dim The dimensionality of the FFT to be transformed"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftSetPlanDim(plHandle: clfftPlanHandle, dim: clfftDim) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the length of each dimension of the FFT"]
    #[doc = "  @details The user must pass a reference to a size_t array, which is set to the"]
    #[doc = "  length of each discrete dimension of the FFT"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dim Dimension of the FFT; describes how many elements are in the clLengths array"]
    #[doc = "  @param[out] clLengths An array of length of size 'dim';  each array value describes the length of each dimension"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftGetPlanLength(
        plHandle: clfftPlanHandle,
        dim: clfftDim,
        clLengths: *mut usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the length of each dimension of the FFT"]
    #[doc = "  @details Sets the plan property which is the length of each discrete dimension of the FFT"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dim The dimension of the FFT; describes how many elements are in the clLengths array"]
    #[doc = "  @param[in] clLengths An array of length of size 'dim';  each array value describes the length of each dimension"]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftSetPlanLength(
        plHandle: clfftPlanHandle,
        dim: clfftDim,
        clLengths: *const usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the distance between consecutive elements of input buffers in each dimension."]
    #[doc = "  @details Depending on how the dimension is set in the plan (for 2D or 3D FFT), strideY or strideZ can be safely"]
    #[doc = "  ignored"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dim The dimension of the stride parameters; provides the number of elements in the array"]
    #[doc = "  @param[out] clStrides An array of strides, of size 'dim'."]
    pub fn clfftGetPlanInStride(
        plHandle: clfftPlanHandle,
        dim: clfftDim,
        clStrides: *mut usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the distance between consecutive elements of input buffers in each dimension."]
    #[doc = "  @details Set the plan properties which is the distance between elements in all dimensions of the input buffer"]
    #[doc = "  (units are in terms of clfftPrecision)"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dim The dimension of the stride parameters; provides the number of elements in the clStrides array"]
    #[doc = "  @param[in] clStrides An array of strides of size 'dim'. Usually, strideX=1 so that successive elements in the first dimension are stored contiguously."]
    #[doc = " \tTypically, strideY=LenX and strideZ=LenX*LenY with the successive elements in the second and third dimensions stored in packed format."]
    #[doc = "  See  @ref DistanceStridesandPitches for details."]
    pub fn clfftSetPlanInStride(
        plHandle: clfftPlanHandle,
        dim: clfftDim,
        clStrides: *mut usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the distance between consecutive elements of output buffers in each dimension."]
    #[doc = "  @details Depending on how the dimension is set in the plan (for 2D or 3D FFT), strideY or strideZ can be safely"]
    #[doc = "  ignored"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dim The dimension of the stride parameters; provides the number of elements in the clStrides array"]
    #[doc = "  @param[out] clStrides An array of strides, of size 'dim'."]
    pub fn clfftGetPlanOutStride(
        plHandle: clfftPlanHandle,
        dim: clfftDim,
        clStrides: *mut usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the distance between consecutive elements of output buffers in a dimension."]
    #[doc = "  @details Sets the plan properties which is the distance between elements in all dimensions of the output buffer"]
    #[doc = "  (units are in terms of clfftPrecision)"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dim The dimension of the stride parameters; provides the number of elements in the clStrides array"]
    #[doc = "  @param[in] clStrides An array of strides of size 'dim'.  Usually, strideX=1 so that successive elements in the first dimension are stored contiguously."]
    #[doc = " \tTypically, strideY=LenX and strideZ=LenX*LenY cause the successive elements in the second and third dimensions be stored in packed format."]
    #[doc = "  @sa clfftSetPlanInStride"]
    pub fn clfftSetPlanOutStride(
        plHandle: clfftPlanHandle,
        dim: clfftDim,
        clStrides: *mut usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the distance between array objects"]
    #[doc = "  @details Pitch is the distance between each discrete array object in an FFT array. This is only used"]
    #[doc = "  for 'array' dimensions in clfftDim; see clfftSetPlanDimension (units are in terms of clfftPrecision)"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] iDist The distance between the beginning elements of the discrete array objects in input buffer."]
    #[doc = "  For contiguous arrays in memory, iDist=(strideX*strideY*strideZ)"]
    #[doc = "  @param[out] oDist The distance between the beginning elements of the discrete array objects in output buffer."]
    #[doc = "  For contiguous arrays in memory, oDist=(strideX*strideY*strideZ)"]
    pub fn clfftGetPlanDistance(
        plHandle: clfftPlanHandle,
        iDist: *mut usize,
        oDist: *mut usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the distance between array objects"]
    #[doc = "  @details Pitch is the distance between each discrete array object in an FFT array. This is only used"]
    #[doc = "  for 'array' dimensions in clfftDim; see clfftSetPlanDimension (units are in terms of clfftPrecision)"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] iDist The distance between the beginning elements of the discrete array objects in input buffer."]
    #[doc = "  For contiguous arrays in memory, iDist=(strideX*strideY*strideZ)"]
    #[doc = "  @param[out] oDist The distance between the beginning elements of the discrete array objects in output buffer."]
    #[doc = "  For contiguous arrays in memory, oDist=(strideX*strideY*strideZ)"]
    pub fn clfftSetPlanDistance(
        plHandle: clfftPlanHandle,
        iDist: usize,
        oDist: usize,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the expected layout of the input and output buffers"]
    #[doc = "  @details Input and output buffers can be filled with either Hermitian, complex, or real numbers.  Complex numbers are stored"]
    #[doc = "  in various layouts; this function retrieves the layouts used by input and output"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] iLayout Indicates how the input buffers are laid out in memory"]
    #[doc = "  @param[out] oLayout Indicates how the output buffers are laid out in memory"]
    pub fn clfftGetLayout(
        plHandle: clfftPlanHandle,
        iLayout: *mut clfftLayout,
        oLayout: *mut clfftLayout,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the expected layout of the input and output buffers"]
    #[doc = "  @details Input and output buffers can be filled with either Hermitian, complex, or real numbers.  Complex numbers can be stored"]
    #[doc = "  in various layouts; this function informs the library what layouts to use for input and output"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] iLayout Indicates how the input buffers are laid out in memory"]
    #[doc = "  @param[in] oLayout Indicates how the output buffers are laid out in memory"]
    pub fn clfftSetLayout(
        plHandle: clfftPlanHandle,
        iLayout: clfftLayout,
        oLayout: clfftLayout,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve whether the input buffers are to be overwritten with results"]
    #[doc = "  @details If the setting performs an in-place transform, the input buffers are overwritten with the results of the"]
    #[doc = "  transform.  If the setting performs an out-of-place transforms, the library looks for separate output buffers"]
    #[doc = "  on the Enqueue call."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] placeness Informs the library to either overwrite the input buffers with results or to write them in separate output buffers"]
    pub fn clfftGetResultLocation(
        plHandle: clfftPlanHandle,
        placeness: *mut clfftResultLocation,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set whether the input buffers are to be overwritten with results"]
    #[doc = "  @details If the setting performs an in-place transform, the input buffers are overwritten with the results of the"]
    #[doc = "  transform.  If the setting performs an out-of-place transforms, the library looks for separate output buffers"]
    #[doc = "  on the Enqueue call."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] placeness Informs the library to either overwrite the input buffers with results or to write them in separate output buffers"]
    pub fn clfftSetResultLocation(
        plHandle: clfftPlanHandle,
        placeness: clfftResultLocation,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Retrieve the final transpose setting of a multi-dimensional FFT"]
    #[doc = "  @details A multi-dimensional FFT transposes the data several times during calculation. If the client"]
    #[doc = "  does not care about the final transpose, to put data back in proper dimension, the final transpose can be skipped"]
    #[doc = "  to improve speed"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] transposed Specifies whether the final transpose can be skipped"]
    pub fn clfftGetPlanTransposeResult(
        plHandle: clfftPlanHandle,
        transposed: *mut clfftResultTransposed,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Set the final transpose setting of a multi-dimensional FFT"]
    #[doc = "  @details A multi-dimensional FFT transposes the data several times during calculation.  If the client"]
    #[doc = "  does not care about the final transpose, to put data back in proper dimension, the final transpose can be skipped"]
    #[doc = "  to improve speed"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] transposed Specifies whether the final transpose can be skipped"]
    pub fn clfftSetPlanTransposeResult(
        plHandle: clfftPlanHandle,
        transposed: clfftResultTransposed,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Get buffer size (in bytes), which may be needed internally for an intermediate buffer"]
    #[doc = "  @details Very large FFT transforms may need multiple passes, and the operation needs a temporary buffer to hold"]
    #[doc = "  intermediate results. This function is only valid after the plan is baked, otherwise, an invalid operation error"]
    #[doc = "  is returned. If the returned buffersize is 0, the runtime needs no temporary buffer."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[out] buffersize Size in bytes for intermediate buffer"]
    pub fn clfftGetTmpBufSize(plHandle: clfftPlanHandle, buffersize: *mut usize) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Register the callback parameters"]
    #[doc = "  @details Client can provide a callback function to do custom processing while reading input data and/or"]
    #[doc = "  writing output data. The callback function is provided as a string."]
    #[doc = "  clFFT library incorporates the callback function string into the main FFT kernel. This function is used"]
    #[doc = "  by client to set the necessary parameters for callback"]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] funcName Callback function name"]
    #[doc = "  @param[in] funcString Callback function in string form"]
    #[doc = "  @param[in] localMemSize Optional - Size (bytes) of the local memory used by callback function; pass 0 if no local memory is used"]
    #[doc = "  @param[in] callbackType Type of callback - Pre-Callback or Post-Callback"]
    #[doc = "  @param[in] userdata Supplementary data if any used by callback function"]
    #[doc = "  @param[in] numUserdataBuffers Number of userdata buffers"]
    pub fn clfftSetPlanCallback(
        plHandle: clfftPlanHandle,
        funcName: *const ::std::os::raw::c_char,
        funcString: *const ::std::os::raw::c_char,
        localMemSize: ::std::os::raw::c_int,
        callbackType: clfftCallbackType,
        userdata: *mut cl_mem,
        numUserdataBuffers: ::std::os::raw::c_int,
    ) -> clfftStatus;
}
extern "C" {
    #[doc = " @brief Enqueue an FFT transform operation, and return immediately (non-blocking)"]
    #[doc = "  @details This transform API function computes the FFT transform. It is non-blocking as it"]
    #[doc = "  only enqueues the OpenCL kernels for execution. The synchronization step must be managed by the user."]
    #[doc = "  @param[in] plHandle Handle to a previously created plan"]
    #[doc = "  @param[in] dir Forward or backward transform"]
    #[doc = "  @param[in] numQueuesAndEvents Number of command queues in commQueues; number of expected events to be returned in outEvents"]
    #[doc = "  @param[in] commQueues An array of cl_command_queues created by the client; the command queues must be a proper subset of"]
    #[doc = " \tthe devices included in the OpenCL context associated with the plan"]
    #[doc = "  @param[in] numWaitEvents Specify the number of elements in the eventWaitList array"]
    #[doc = "  @param[in] waitEvents Events for which the transform waits to complete before executing on the device"]
    #[doc = "  @param[out] outEvents The runtime fills this array with events corresponding one to one with the input command queues passed"]
    #[doc = "\tin commQueues.  This parameter can have the value NULL or nullptr. When the value is NULL, the client is not interested in receiving notifications"]
    #[doc = "\twhen transforms are finished, otherwise, (if not NULL) the client is responsible for allocating this array with at least"]
    #[doc = "\tas many elements as specified in numQueuesAndEvents."]
    #[doc = "  @param[in] inputBuffers An array of cl_mem objects that contain data for processing by the FFT runtime. If the transform"]
    #[doc = "  is in-place, the FFT results overwrite the input buffers"]
    #[doc = "  @param[out] outputBuffers An array of cl_mem objects that store the results of out-of-place transforms. If the transform"]
    #[doc = "  is in-place, this parameter may be NULL or nullptr and is completely ignored"]
    #[doc = "  @param[in] tmpBuffer A cl_mem object that is reserved as a temporary buffer for FFT processing. If clTmpBuffers is NULL or nullptr,"]
    #[doc = "  and the library needs temporary storage, an internal temporary buffer is created on the fly managed by the library."]
    #[doc = "  @return Enum describing error condition; superset of OpenCL error codes"]
    pub fn clfftEnqueueTransform(
        plHandle: clfftPlanHandle,
        dir: clfftDirection,
        numQueuesAndEvents: cl_uint,
        commQueues: *mut cl_command_queue,
        numWaitEvents: cl_uint,
        waitEvents: *const cl_event,
        outEvents: *mut cl_event,
        inputBuffers: *mut cl_mem,
        outputBuffers: *mut cl_mem,
        tmpBuffer: cl_mem,
    ) -> clfftStatus;
}
