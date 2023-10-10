use std::io::{Write, Read, Error, ErrorKind};

pub struct Encoder8 {}

impl Encoder8 {
    pub fn encode<Sink>(x: u64, sink: &mut Sink) -> std::io::Result<()>
    where Sink: Write {
        if x <= 127 {
            return sink.write_all(&[(x as u8) << 1]);
        }

        let num_bits: u32 = 64 - x.leading_zeros();
        let num_bytes: usize = ((num_bits + 6) / 7) as usize;
        let mut buf: [u8; 10] = [0; 10];
        let mut y = x;
        for i in 0..num_bytes {
            buf[i] = ((y & 0x7f) as u8) << 1;
            y >>= 7;
        }
        let mut res: [u8; 10] = [0; 10];
        let mut i = num_bytes - 1;
        let mut j = 0;
        while i > 0 {
            res[j] = buf[i] | 1;
            i -= 1;
            j += 1;
        }
        res[j] = buf[0];
        sink.write_all(&res[0..num_bytes])
    }
}

pub struct Decoder8 {}

impl Decoder8 {
    pub fn decode<'a, Source>(source: &mut Source) -> std::io::Result<Option<u64>>
    where
        Source: Read,
    {
        let mut buf: [u8; 1] = [0; 1];
        let mut count = source.read(&mut buf)?;
        match count {
            0 => Ok(None),
            _ => {
                let mut b = buf[0];
                let mut x = (b as u64) >> 1;
                let mut v = b;
                while v & 1 == 1 {
                    count = source.read(&mut buf)?;
                    if count == 1 {
                        b = buf[0];
                        x = (x << 7) | ((b >> 1) as u64);
                        v = b;
                    } else {
                        return Err(Error::from(ErrorKind::UnexpectedEof));
                    }
                }
                Ok(Some(x))
            }
        }
    }
}

#[allow(dead_code)]
pub struct Encoder16 {}

impl Encoder16 {
    #[allow(dead_code)]
    pub fn encode(x: u64, sink: &mut Vec<u16>) {
        if x <= 0x7fff {
            sink.push((x << 1) as u16);
            return;
        }

        let num_bits: u32 = 64 - x.leading_zeros();
        let num_bytes: usize = ((num_bits + 14) / 15) as usize;
        let mut buf: [u16; 5] = [0, 0, 0, 0, 0];
        let mut y = x;
        for i in 0..num_bytes {
            buf[i] = ((y & 0x7fff) as u16) << 1;
            y >>= 15;
        }
        let mut i = num_bytes - 1;
        while i > 0 {
            sink.push(buf[i] | 1);
            i -= 1;
        }
        sink.push(buf[0]);
    }
}

#[allow(dead_code)]
pub struct Decoder16 {}

impl Decoder16 {
    #[allow(dead_code)]
    pub fn decode<'a, I>(iter: &mut I) -> Option<u64>
    where
        I: Iterator<Item = &'a u16>,
    {
        match iter.next() {
            None => None,
            Some(b) => {
                let mut x = (*b as u64) >> 1;
                let mut v = *b;
                while v & 1 == 1 {
                    if let Some(w) = iter.next() {
                        x = (x << 15) | ((*w >> 1) as u64);
                        v = *w;
                    } else {
                        return None;
                    }
                }
                Some(x)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn vbyte_test_1() {
        let x: u64 = 0;
        let mut bytes: Vec<u8> = Vec::new();
        Encoder8::encode(x, &mut bytes).expect("encode failed");
        assert_eq!(bytes, vec![0]);
        let mut source = Cursor::new(bytes);
        let res = Decoder8::decode(&mut source).expect("decode failed");
        assert_eq!(res, Some(x));
    }

    #[test]
    fn vbyte_test_2() {
        let x: u64 = 127;
        let mut bytes: Vec<u8> = Vec::new();
        Encoder8::encode(x, &mut bytes).expect("encode failed");
        assert_eq!(bytes, vec![254]);
        let mut source = Cursor::new(bytes);
        let res = Decoder8::decode(&mut source).expect("decode failed");
        assert_eq!(res, Some(x));
    }

    #[test]
    fn vbyte_test_3() {
        let x: u64 = 128;
        let mut bytes: Vec<u8> = Vec::new();
        Encoder8::encode(x, &mut bytes).expect("encode failed");
        assert_eq!(bytes, vec![3, 0]);
        let mut source = Cursor::new(bytes);
        let res = Decoder8::decode(&mut source).expect("decode failed");
        assert_eq!(res, Some(x));
    }

    #[test]
    fn vbyte_test_4() {
        let xs: [u64; 16] = [
            0xf,
            0x54,
            0xb99,
            0xe6f9,
            0x33fa1,
            0x212ce9,
            0x7e496e2,
            0xaae84a39,
            0x5f89b315d,
            0xa8b333445b,
            0x2179bec851b,
            0xcf21a080f954,
            0x82a42d77441ec,
            0x31c12a609978ed,
            0x3f7ce900e5bfc6f,
            0x6de3d6138095fdd9,
        ];
        for x in xs {
            let mut bytes: Vec<u8> = Vec::new();
            Encoder8::encode(x, &mut bytes).expect("encode failed");
            let mut source = Cursor::new(bytes);
            let res = Decoder8::decode(&mut source).expect("decode failed");
            assert_eq!(res, Some(x));
        }
    }
}
