#include <iostream>
#include <stdio.h>
#include <sys/time.h>  /* gettimeofday */
// #include <assert.h>
// #include <string.h> /* memcpy */
#include <jpeglib.h>

#if 0
typedef enum FilterMode {
  kFilterNone = 0,      // Point sample; Fastest.
  kFilterLinear = 1,    // Filter horizontally only.
  kFilterBilinear = 2,  // Faster than box, but lower quality scaling down.
  kFilterBox = 3        // Highest quality.
} FilterModeEnum;

static __inline int Abs(int v) {
  return v >= 0 ? v : -v;
}

// Blend 2 rows into 1.
static inline void HalfRow_C(const uint8_t* src_uv,
                      ptrdiff_t src_uv_stride,
                      uint8_t* dst_uv,
                      int width) {
  int x;
  for (x = 0; x < width; ++x) {
    dst_uv[x] = (src_uv[x] + src_uv[src_uv_stride + x] + 1) >> 1;
  }
}

// C version 2x2 -> 2x1.
void InterpolateRow_C(uint8_t* dst_ptr,
                      const uint8_t* src_ptr,
                      ptrdiff_t src_stride,
                      int width,
                      int source_y_fraction) {
  int y1_fraction = source_y_fraction;
  int y0_fraction = 256 - y1_fraction;
  const uint8_t* src_ptr1 = src_ptr + src_stride;
  int x;
  assert(source_y_fraction >= 0);
  assert(source_y_fraction < 256);

  if (y1_fraction == 0) {
    memcpy(dst_ptr, src_ptr, width);
    return;
  }
  if (y1_fraction == 128) {
    HalfRow_C(src_ptr, src_stride, dst_ptr, width);
    return;
  }
  for (x = 0; x < width; ++x) {
    dst_ptr[0] =
        (src_ptr[0] * y0_fraction + src_ptr1[0] * y1_fraction + 128) >> 8;
    ++src_ptr;
    ++src_ptr1;
    ++dst_ptr;
  }
}

#define BLENDER1(a, b, f) ((a) * (0x7f ^ f) + (b)*f) >> 7
#define BLENDERC(a, b, f, s) \
  (uint32_t)(BLENDER1(((a) >> s) & 255, ((b) >> s) & 255, f) << s)
#define BLENDER(a, b, f)                                                 \
  BLENDERC(a, b, f, 24) | BLENDERC(a, b, f, 16) | BLENDERC(a, b, f, 8) | \
      BLENDERC(a, b, f, 0)

void ScaleFilterCols_C(uint8_t* dst_ptr,
                       const uint8_t* src_ptr,
                       int dst_width,
                       int x,
                       int dx) {
  int j;
  for (j = 0; j < dst_width - 1; j += 2) {
    int xi = x >> 16;
    int a = src_ptr[xi];
    int b = src_ptr[xi + 1];
    dst_ptr[0] = BLENDER(a, b, x & 0xffff);
    x += dx;
    xi = x >> 16;
    a = src_ptr[xi];
    b = src_ptr[xi + 1];
    dst_ptr[1] = BLENDER(a, b, x & 0xffff);
    x += dx;
    dst_ptr += 2;
  }
  if (dst_width & 1) {
    int xi = x >> 16;
    int a = src_ptr[xi];
    int b = src_ptr[xi + 1];
    dst_ptr[0] = BLENDER(a, b, x & 0xffff);
  }
}

#define align_buffer_64(var, size)                                           \
  uint8_t* var##_mem = (uint8_t*)(malloc((size) + 63));         /* NOLINT */ \
  uint8_t* var = (uint8_t*)(((intptr_t)(var##_mem) + 63) & ~63) /* NOLINT */

#define free_aligned_buffer_64(var) \
  free(var##_mem);                  \
  var = 0

// Divide num by div and return as 16.16 fixed point result.
int FixedDiv(int num, int div) {
  return (int)(((int64_t)(num) << 16) / div);
}

// Divide num - 1 by div - 1 and return as 16.16 fixed point result.
int FixedDiv1(int num, int div) {
  return (int)((((int64_t)(num) << 16) - 0x00010001) / (div - 1));
}

#define CENTERSTART(dx, s) (dx < 0) ? -((-dx >> 1) + s) : ((dx >> 1) + s)

// Compute slope values for stepping.
void ScaleSlope(int src_width,
                int src_height,
                int dst_width,
                int dst_height,
                enum FilterMode filtering,
                int* x,
                int* y,
                int* dx,
                int* dy) {
  assert(x != NULL);
  assert(y != NULL);
  assert(dx != NULL);
  assert(dy != NULL);
  assert(src_width != 0);
  assert(src_height != 0);
  assert(dst_width > 0);
  assert(dst_height > 0);
  // Check for 1 pixel and avoid FixedDiv overflow.
  if (dst_width == 1 && src_width >= 32768) {
    dst_width = src_width;
  }
  if (dst_height == 1 && src_height >= 32768) {
    dst_height = src_height;
  }
  if (filtering == kFilterBox) {
    // Scale step for point sampling duplicates all pixels equally.
    *dx = FixedDiv(Abs(src_width), dst_width);
    *dy = FixedDiv(src_height, dst_height);
    *x = 0;
    *y = 0;
  } else if (filtering == kFilterBilinear) {
    // Scale step for bilinear sampling renders last pixel once for upsample.
    if (dst_width <= Abs(src_width)) {
      *dx = FixedDiv(Abs(src_width), dst_width);
      *x = CENTERSTART(*dx, -32768);  // Subtract 0.5 (32768) to center filter.
    } else if (src_width > 1 && dst_width > 1) {
      *dx = FixedDiv1(Abs(src_width), dst_width);
      *x = 0;
    }
    if (dst_height <= src_height) {
      *dy = FixedDiv(src_height, dst_height);
      *y = CENTERSTART(*dy, -32768);  // Subtract 0.5 (32768) to center filter.
    } else if (src_height > 1 && dst_height > 1) {
      *dy = FixedDiv1(src_height, dst_height);
      *y = 0;
    }
  } else if (filtering == kFilterLinear) {
    // Scale step for bilinear sampling renders last pixel once for upsample.
    if (dst_width <= Abs(src_width)) {
      *dx = FixedDiv(Abs(src_width), dst_width);
      *x = CENTERSTART(*dx, -32768);  // Subtract 0.5 (32768) to center filter.
    } else if (src_width > 1 && dst_width > 1) {
      *dx = FixedDiv1(Abs(src_width), dst_width);
      *x = 0;
    }
    *dy = FixedDiv(src_height, dst_height);
    *y = *dy >> 1;
  } else {
    // Scale step for point sampling duplicates all pixels equally.
    *dx = FixedDiv(Abs(src_width), dst_width);
    *dy = FixedDiv(src_height, dst_height);
    *x = CENTERSTART(*dx, 0);
    *y = CENTERSTART(*dy, 0);
  }
  // Negative src_width means horizontally mirror.
  if (src_width < 0) {
    *x += (dst_width - 1) * *dx;
    *dx = -*dx;
    // src_width = -src_width;   // Caller must do this.
  }
}
#undef CENTERSTART

void ScalePlaneBilinearUp(int src_width,
                          int src_height,
                          int dst_width,
                          int dst_height,
                          int src_stride,
                          int dst_stride,
                          const uint8_t* src_ptr,
                          uint8_t* dst_ptr,
                          enum FilterMode filtering) {
  int j;
  // Initial source x/y coordinate and step values as 16.16 fixed point.
  int x = 0;
  int y = 0;
  int dx = 0;
  int dy = 0;
  const int max_y = (src_height - 1) << 16;
  void (*InterpolateRow)(uint8_t * dst_ptr, const uint8_t* src_ptr,
                         ptrdiff_t src_stride, int dst_width,
                         int source_y_fraction) = InterpolateRow_C;
  void (*ScaleFilterCols)(uint8_t * dst_ptr, const uint8_t* src_ptr,
                          int dst_width, int x, int dx) = ScaleFilterCols_C;

  ScaleSlope(src_width, src_height, dst_width, dst_height, filtering, &x, &y,
             &dx, &dy);
  src_width = Abs(src_width);

  if (y > max_y) {
    y = max_y;
  }

  {
    int yi = y >> 16;
    const uint8_t* src = src_ptr + yi * (int64_t)src_stride;

    // Allocate 2 row buffers.
    const int kRowSize = (dst_width + 31) & ~31;
    align_buffer_64(row, kRowSize * 2);

    uint8_t* rowptr = row;
    int rowstride = kRowSize;
    int lasty = yi;

    ScaleFilterCols(rowptr, src, dst_width, x, dx);
    if (src_height > 1) {
      src += src_stride;
    }
    ScaleFilterCols(rowptr + rowstride, src, dst_width, x, dx);
    if (src_height > 2) {
      src += src_stride;
    }

    for (j = 0; j < dst_height; ++j) {
      yi = y >> 16;
      if (yi != lasty) {
        if (y > max_y) {
          y = max_y;
          yi = y >> 16;
          src = src_ptr + yi * (int64_t)src_stride;
        }
        if (yi != lasty) {
          ScaleFilterCols(rowptr, src, dst_width, x, dx);
          rowptr += rowstride;
          rowstride = -rowstride;
          lasty = yi;
          if ((y + 65536) < max_y) {
            src += src_stride;
          }
        }
      }

      if (filtering == kFilterLinear) {
        InterpolateRow(dst_ptr, rowptr, 0, dst_width, 0);
      } else {
        int yf = (y >> 8) & 255;
        InterpolateRow(dst_ptr, rowptr, rowstride, dst_width, yf);
      }
      dst_ptr += dst_stride;
      y += dy;
    }
    free_aligned_buffer_64(row);
  }
}
#endif

void resize_plane(uint8_t *src, uint8_t *dst, int srcWidth, int srcHeight, int dstWidth, int dstHeight)
{
	int sw = srcWidth;
	int sh = srcHeight;
	int dw = dstWidth;
	int dh = dstHeight;
	int y, x;
	unsigned long int srcy, srcx;
	unsigned long int xrIntFloat_16 = (sw << 16) / dw + 1;
	unsigned long int yrIntFloat_16 = (sh << 16) / dh + 1;

	uint8_t *dst_y_slice = dst;
	uint8_t *src_y_slice;

	for (y = 0; y < (dh & ~1); ++y)
	{
		srcy = (y * yrIntFloat_16) >> 16;
		src_y_slice = src + srcy * sw;

		for(x = 0; x < (dw & ~1); ++x)
		{
			srcx = (x * xrIntFloat_16) >> 16;
			dst_y_slice[x] = src_y_slice[srcx];
		}
		dst_y_slice += dw;
	}
}

void resize_uv_plane(uint8_t *src_u, uint8_t *src_v, uint8_t *dst, int srcWidth, int srcHeight, int dstWidth, int dstHeight)
{
	int sw = srcWidth;
	int sh = srcHeight;
	int dw = dstWidth;
	int dh = dstHeight;
	int y, x;
	unsigned long int srcy, srcx;
	unsigned long int xrIntFloat_16 = (sw << 16) / dw + 1;
	unsigned long int yrIntFloat_16 = (sh << 16) / dh + 1;

	uint8_t *dst_uv_slice = dst;
	uint8_t *src_u_slice = src_u;
	uint8_t *src_v_slice = src_v;

	for (y = 0; y < (dh & ~1); ++y)
	{
		srcy = (y * yrIntFloat_16) >> 16;
		src_u_slice = src_u + srcy * sw;
		src_v_slice = src_v + srcy * sw;

		for(x = 0; x < (dw & ~1); ++x)
		{
			srcx = (x * xrIntFloat_16) >> 16;
			dst_uv_slice[2*x] 	= src_u_slice[srcx];
			dst_uv_slice[2*x+1] = src_v_slice[srcx];
		}
		dst_uv_slice += (dw*2);
	}
}

#define src_width 2560
#define src_height 1440
/* 缩放2.5倍 */
// #define resize_width 6400
// #define resize_height 3600
// #define resize_height_pre_proc 40

/* 缩放2倍 */
// #define resize_width 5120
// #define resize_height 2880
// #define resize_height_pre_proc 32

/* 缩放1.5倍 */
// #define resize_width 3840
// #define resize_height 2160
// #define resize_height_pre_proc 24

/* 宽 缩放2.3倍 高缩放3倍 */
// #define resize_width 5888
// #define resize_height 4320
// #define resize_height_pre_proc 48

/* 宽 缩放2.8倍 高缩放3倍 整体3000W */
#define resize_width 7168
#define resize_height 4320
#define resize_height_pre_proc 48

// #define SAVE_JPEG_TO_YUV_SRC_FILE
// #define TEST_SCALE_JPEG

void jpegtoYUV()
{
  struct jpeg_decompress_struct info_;
  struct jpeg_error_mgr e_;
  info_.err = jpeg_std_error(&e_);
  jpeg_create_decompress(&info_);

  FILE *infile = fopen("test.jpg", "rb");
  jpeg_stdio_src(&info_, infile);

#ifdef TEST_SCALE_JPEG
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  FILE *outfile = fopen("./jpeg_scale_end.jpg", "wb");
  jpeg_stdio_dest(&cinfo, outfile);
#endif

  // Get all compression parameters.
  jpeg_read_header(&info_, 1);

  // Configure the decompressors for raw data output.
  info_.out_color_space = JCS_YCbCr;
  // info_.do_fancy_upsampling = TRUE;
  info_.raw_data_out = TRUE;
  jpeg_start_decompress(&info_);

#ifdef TEST_SCALE_JPEG
  cinfo.image_width = info_.output_width * 2;
  cinfo.image_height = info_.output_height * 2;
  cinfo.input_components = info_.output_components;

  jpeg_set_defaults(&cinfo);
  // cinfo.optimize_coding = TRUE;
  cinfo.dct_method = JDCT_IFAST;
  jpeg_set_quality(&cinfo, 50, TRUE /* limit to baseline-JPEG values */);

  cinfo.raw_data_in = TRUE;
  cinfo.in_color_space = info_.out_color_space;
  cinfo.comp_info[0].h_samp_factor = 2;
  cinfo.comp_info[0].v_samp_factor = 2;

  jpeg_start_compress(&cinfo, TRUE);
#endif

  int nPicHeight = info_.output_height;
  int nPicWidth  = info_.output_width;

#ifdef SAVE_JPEG_TO_YUV_SRC_FILE
  int nLen = nPicHeight * nPicWidth * 3 / 2;
  unsigned char *yuv420p_buf = (unsigned char *)malloc(nLen);
  unsigned char *m_pYbuffer = yuv420p_buf;
  unsigned char *m_pUbuffer = yuv420p_buf + nPicHeight * nPicWidth;
  unsigned char *m_pVbuffer = m_pUbuffer + (nPicHeight * nPicWidth >> 2);
#else
  int nLen = nPicWidth * 24;
  unsigned char *yuv420p_buf = (unsigned char *)malloc(nLen);
  unsigned char *m_pYbuffer = yuv420p_buf;
  unsigned char *m_pUbuffer = yuv420p_buf + (nPicWidth << 4);
  unsigned char *m_pVbuffer = m_pUbuffer + (nPicWidth << 2);
#endif

#ifdef TEST_SCALE_JPEG
  unsigned char *yuv420p_buf_scaled = (unsigned char *)malloc(nPicWidth * 32 * 3);
  unsigned char *m_pYbuffer_scaled = yuv420p_buf_scaled;
  unsigned char *m_pUbuffer_scaled = m_pYbuffer_scaled + (nPicWidth << 6);
  unsigned char *m_pVbuffer_scaled = m_pUbuffer_scaled + (nPicWidth << 4);
#else
  unsigned char *yuv420p_buf_scaled = (unsigned char *)malloc(resize_width * resize_height * 3 / 2);
  unsigned char *m_pYbuffer_scaled = yuv420p_buf_scaled;
  unsigned char *m_pUbuffer_scaled = m_pYbuffer_scaled + resize_width * resize_height;
  unsigned char *m_pVbuffer_scaled = m_pUbuffer_scaled + resize_width * resize_height / 4;
#endif

  //分配内存获取数据
  int ci;
  // int yx, yr,ci,bi;
  JSAMPIMAGE image_one = (JSAMPIMAGE)malloc(3 * sizeof(JSAMPARRAY));

#ifdef TEST_SCALE_JPEG
  JSAMPIMAGE image_dst_0 = (JSAMPIMAGE)malloc(3 * sizeof(JSAMPARRAY));
  JSAMPIMAGE image_dst_1 = (JSAMPIMAGE)malloc(3 * sizeof(JSAMPARRAY));
#endif

  for (ci = 0; ci < 3; ci++)
  {
      image_one[ci] = (JSAMPARRAY)malloc(16 * sizeof(JSAMPROW));
#ifdef TEST_SCALE_JPEG
      image_dst_0[ci] = (JSAMPARRAY)malloc(16 * sizeof(JSAMPROW));
      image_dst_1[ci] = (JSAMPARRAY)malloc(16 * sizeof(JSAMPROW));
#endif
  }

  unsigned char *ytmp_src = m_pYbuffer;
  unsigned char *utmp_src = m_pUbuffer;
  unsigned char *vtmp_src = m_pVbuffer;

  /* 读取完整YUV数据，减少内存拷贝 */
  for (int i = 0; i < 16; ++i, ytmp_src += info_.output_width)
  {
    image_one[0][i] = ytmp_src; //y 分量空间初始化
  }

  for (int i = 0; i < 16; i += 2, utmp_src += info_.output_width / 2, vtmp_src += info_.output_width / 2)
  {
    image_one[1][i / 2] = utmp_src; //u 分量初始化
    image_one[2][i / 2] = vtmp_src; //v 分量初始化
  }

  unsigned char *ytmp = m_pYbuffer_scaled;
  unsigned char *utmp = m_pUbuffer_scaled;
  unsigned char *vtmp = m_pVbuffer_scaled;

#ifdef TEST_SCALE_JPEG
  /* 分段读取数据、缩放、编码 */
  for (int i = 0; i < 16; ++i, ytmp += cinfo.image_width)
  {
    image_dst_0[0][i] = ytmp;
    image_dst_1[0][i] = ytmp + 16 * cinfo.image_width;
  }

  for (int i = 0; i < 16; i += 2, utmp += nPicWidth, vtmp += nPicWidth)
  {
    image_dst_0[1][i / 2] = utmp;
    image_dst_1[1][i / 2] = utmp + 8 * nPicWidth;

    image_dst_0[2][i / 2] = vtmp;
    image_dst_1[2][i / 2] = vtmp + 8 * nPicWidth;
  }
#endif

  for (int i = 0; i < info_.output_height; i += 16)
  {
    int nRows = 16;
    // if ((info_.output_height) < (i + 16))
    // {
    //     nRows = info_.output_height - i;
    // }

    jpeg_read_raw_data(&info_, image_one, (JDIMENSION)nRows);

    /* 按通道缩放 */
    // ScalePlaneBilinearUp(nPicWidth, nRows,  nPicWidth * 2, nRows * 2, nPicWidth, nPicWidth*2, m_pYbuffer, ytmp, kFilterBilinear);
    // ScalePlaneBilinearUp(nPicWidth, nRows/4, nPicWidth * 2, nRows / 2, nPicWidth, nPicWidth*2, m_pUbuffer, utmp, kFilterBilinear);
    // ScalePlaneBilinearUp(nPicWidth, nRows/4, nPicWidth * 2, nRows / 2, nPicWidth, nPicWidth*2, m_pVbuffer, vtmp, kFilterBilinear);

#ifdef TEST_SCALE_JPEG
    resize_plane(m_pYbuffer, m_pYbuffer_scaled, nPicWidth, nRows,  nPicWidth * 2, nRows * 2);
    resize_plane(m_pUbuffer, m_pUbuffer_scaled, nPicWidth/2, nRows/2,  nPicWidth, nRows);
    resize_plane(m_pVbuffer, m_pVbuffer_scaled, nPicWidth/2, nRows/2,  nPicWidth, nRows);

    /* 编码为jpeg图片 */
    jpeg_write_raw_data(&cinfo, image_dst_0, nRows);
    jpeg_write_raw_data(&cinfo, image_dst_1, nRows);
#else
    /* save yuv */
    resize_plane(image_one[0][0], ytmp, nPicWidth, nRows,  resize_width,  resize_height_pre_proc);
    // resize_plane(image_one[1][0], utmp, nPicWidth/2, nRows/2,  nPicWidth, nRows);
    // resize_plane(image_one[2][0], vtmp, nPicWidth/2, nRows/2,  nPicWidth, nRows);

    /* save nv21 */
    resize_uv_plane(image_one[1][0], image_one[2][0], utmp, nPicWidth/2, nRows/2,  resize_width / 2, resize_height_pre_proc / 2);

    /* 目的YUV数据偏移(整张图片) */
    ytmp += (resize_width * resize_height_pre_proc);
    /*save yuv*/
    // utmp += (nPicWidth * nRows);
    // vtmp += (nPicWidth * nRows);

    /* save nv21 */
    utmp += (resize_width * resize_height_pre_proc / 2);
#endif

#ifdef SAVE_JPEG_TO_YUV_SRC_FILE
  /* 读取完整YUV数据，减少内存拷贝 */
  for (int i = 0; i < 16; ++i, ytmp_src += info_.output_width)
  {
    image_one[0][i] = ytmp_src; //y 分量空间初始化
  }

  for (int i = 0; i < 16; i += 2, utmp_src += info_.output_width / 2, vtmp_src += info_.output_width / 2)
  {
    image_one[1][i / 2] = utmp_src; //u 分量初始化
    image_one[2][i / 2] = vtmp_src; //v 分量初始化
  }
#endif
  }

#ifndef TEST_SCALE_JPEG
  FILE *yuv_outfile = fopen("dst.yuv", "wb");
  fwrite(yuv420p_buf_scaled, resize_width * resize_height * 3 / 2, 1, yuv_outfile);
  fflush(yuv_outfile);
  fclose(yuv_outfile);
#endif

#ifdef SAVE_JPEG_TO_YUV_SRC_FILE
  FILE *yuv_srcfile = fopen("jpg_to_yuv.yuv", "wb");
  fwrite(yuv420p_buf, nLen, 1, yuv_outfile);
  fflush(yuv_srcfile);
  fclose(yuv_srcfile);
#endif

  for (ci = 0; ci < 3; ci++)
  {
    free(image_one[ci]);
#ifdef TEST_SCALE_JPEG
    free(image_dst_0[ci]);
    free(image_dst_1[ci]);
#endif
  }
  free(image_one);
#ifdef TEST_SCALE_JPEG
  free(image_dst_0);
  free(image_dst_1);
#endif
  free(yuv420p_buf);
  free(yuv420p_buf_scaled);

  jpeg_finish_decompress(&info_);
  fclose(infile);
  jpeg_destroy_decompress(&info_);

#ifdef TEST_SCALE_JPEG
  jpeg_finish_compress(&cinfo);
  fclose(outfile);
  jpeg_destroy_compress(&cinfo);
#endif

  return;
}

int main()
{
  struct  timeval    ts, te;
  gettimeofday(&ts, NULL);
  std::cout<<ts.tv_sec*1000 + ts.tv_usec/1000 << std::endl;

  jpegtoYUV();
	// test_jpeg_scale();
    //test_libjpeg_scale();
  gettimeofday(&te, NULL);
  std::cout<<te.tv_sec*1000 + te.tv_usec/1000 << std::endl;

  return 0;
}
