use crate::params;

pub struct CloudKey {
    pub decomposition_offset: u32,
}

impl CloudKey {
    pub fn new() -> Self {
        return CloudKey {
            decomposition_offset: gen_decomposition_offset(),
        };
    }
}

pub fn gen_decomposition_offset() -> u32 {
    let mut offset: u32 = 0;

    for i in 0..(params::trgsw_lv1::L as u32) {
        offset = offset.wrapping_add(
            params::trgsw_lv1::BG / 2 * (1 << (32 - (i + 1) * params::trgsw_lv1::BGBIT)),
        );
    }

    return offset;
}
