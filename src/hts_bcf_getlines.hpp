#pragma once
#include "bcfcpp.hpp"

#include <range/v3/detail/prologue.hpp>
namespace ranges
{
    /// \addtogroup group-views
    /// @{
  template<int Fields=BCF_UN_ALL>
    struct getlines_hts_view : view_facade<getlines_hts_view<Fields>, unknown>
    {
    private:
      friend range_access;
      BCFFile * sin_;
      UnpackedBCFLine<Fields> str_;
      struct cursor
      {
      private:
        friend range_access;
        using single_pass = std::true_type;
        getlines_hts_view * rng_ = nullptr;

      public:
        cursor() = default;
        explicit cursor(getlines_hts_view * rng)
          : rng_(rng)
        {};
        void next()
        {
          rng_->next();
        }
        UnpackedBCFLine<Fields> & read() const noexcept
        {
          return rng_->str_;
        }
        bool equal(default_sentinel_t) const
        {
          return !rng_->sin_;
        }
        bool equal(cursor that) const
        {
          return !rng_->sin_ == !that.rng_->sin_;
        }
      };
      void next()
      {
        if(getline_bcf(*sin_,str_)==-1){
          sin_ = nullptr;
        }else{

        }
      }
      cursor begin_cursor()
      {
        return cursor{this};
      }

    public:
      getlines_hts_view() = default;
      getlines_hts_view(BCFFile & sin)
        : sin_(&sin)
        , str_{}
      {
        this->next(); // prime the pump
      }
      BCFLine & cached() noexcept
      {
        return str_;
      }
    };

template <int Fields=BCF_UN_ALL>
  struct getlines_hts_fn
  {
    getlines_hts_view<Fields> operator()(BCFFile & sin, int which=BCF_UN_STR) const
    {
      return getlines_hts_view<Fields>{sin, which};
    }
  };

  RANGES_INLINE_VARIABLE(getlines_hts_fn, getlines_hts)
  /// @}
} // namespace ranges

#include <range/v3/detail/epilogue.hpp>
