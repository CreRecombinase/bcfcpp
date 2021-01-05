#pragma once

#include "htslib/vcf.h"
#include <range/v3/view/facade.hpp>
#include <cstring>
#define bcf_int64_vector_end (-9223372036854775807LL) /* INT64_MIN + 1 */
#define bcf_int64_missing    (-9223372036854775807LL - 1LL)  /* INT64_MIN */
template <class To, class From>
typename std::enable_if_t<
    sizeof(To) == sizeof(From) &&
    std::is_trivially_copyable_v<From> &&
    std::is_trivially_copyable_v<To>,
    To>
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept
{
    static_assert(std::is_trivially_constructible_v<To>,
        "This implementation additionally requires destination type to be trivially constructible");

    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}


template <typename T>
T convert_le(const uint8_t *buf) {
  T ret;
  std::memcpy(&ret,buf,sizeof(T));
  return ret;
}

template <> int8_t convert_le<int8_t>(const uint8_t *buf) { return le_to_i8(buf); }

template <> int16_t convert_le<int16_t>(const uint8_t *buf) { return le_to_i16(buf); }

template<>
int32_t convert_le<int32_t>(const uint8_t *buf){
  return le_to_i32(buf);
}

template<>
float convert_le<float>(const uint8_t *buf){
  return bit_cast<float>(le_to_u32(buf));
}

template <typename T> struct bcf_vector {
  static const T end;
  static const T missing;
};

template <> struct bcf_vector<int8_t> {
  static const int8_t end=bcf_int8_vector_end;
  static const int8_t missing=bcf_int8_missing;
};

template <> struct bcf_vector<char> {
  static const char end=bcf_int8_vector_end;
  static const char missing=bcf_int8_missing;
};

template <> struct bcf_vector<int16_t> {
  static const int16_t end=bcf_int16_vector_end;
  static const int16_t missing=bcf_int16_missing;
};


template <> struct bcf_vector<int32_t> {
  static const int32_t end=bcf_int32_vector_end;
  static const int32_t missing=bcf_int32_missing;
};

template <> struct bcf_vector<int64_t> {
  static const int64_t end=bcf_int64_vector_end;
  static const int64_t missing=bcf_int64_missing;
};

// template <> struct bcf_vector<float> {
//   constexpr static const float  end= std::numeric_limits<float>::signaling_NaN();
//   static const float missing=bcf_float_missing;
// };

//float bcf_float_missing    = memcpy(reinterpret_cast<float>(0x7F800001u);

// HTSLIB_EXPORT
// uint32_t bcf_float_vector_end = 0x7F800002;
template<typename T>
struct always_false : std::false_type {};

template <typename T>
inline bool is_vector_end(uint8_t const* x){
      static_assert(always_false<T>::value , "You must specialize foo<> for your type");
      return false;
}

template <> inline bool is_vector_end<std::int8_t>(uint8_t const *x) {
  return convert_le<std::int8_t>(x) == bcf_int8_vector_end;
}
// template<> inline bool is_vector_end<char>(uint8_t const *x) {
//   return convert_le<char>(x) == bcf_vector<int8_t>::end;
// }
template <> inline bool is_vector_end<std::int16_t>(uint8_t const *x) {
  return convert_le<std::int16_t>(x) == bcf_vector<int16_t>::end;
}
template <> inline bool is_vector_end<std::int32_t>(uint8_t const *x) {
  return convert_le<std::int32_t>(x) == bcf_vector<int32_t>::end;
}
template <> inline bool is_vector_end<std::int64_t>(uint8_t const *x) {
  return convert_le<std::int64_t>(x) == bcf_vector<int64_t>::end;
}


template<>
inline bool is_vector_end<float>(uint8_t const* x) {
  const auto t = std::numeric_limits<float>::signaling_NaN();
  return std::memcmp((void*)&x,(void*)&t,sizeof(float));
}




template<typename T>
class bcf_fmt_conversion_range
  : public ranges::view_facade<bcf_fmt_conversion_range<T>>
{
  friend ranges::range_access;
  uint8_t const *sz_;// = bcf_vector<T>::end;
  T val;
  T const & read() const { return val; }
  bool equal(ranges::default_sentinel_t) const { return is_vector_end<T>(sz_);}
  void next() {
    sz_+=sizeof(T);
    val=convert_le<T>(sz_);
  }
public:
  bcf_fmt_conversion_range() = default;
  explicit bcf_fmt_conversion_range(uint8_t const *sz) : sz_(sz),val(convert_le<T>(sz_))
  {
    assert(sz != nullptr);
  }
};
