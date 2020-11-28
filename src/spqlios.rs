use std::os::raw::c_int;
use std::os::raw::c_uint;
use std::os::raw::c_double;

pub enum SpqliosImpl {}

extern {
    pub fn Spqlios_new(N: c_int) -> *mut SpqliosImpl;
    pub fn Spqlios_destructor(spqlios: *mut SpqliosImpl);
    pub fn Spqlios_ifft_lv1(spqlios: *mut SpqliosImpl, res: *mut c_double, src: *const c_uint);
    pub fn Spqlios_fft_lv1(spqlios: *mut SpqliosImpl, res: *mut c_uint, src: *const c_double);
}