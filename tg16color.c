// MIT License

// Copyright (c) 2020 Jamie Blanks aka Kitrinx

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#define LODEPNG_COMPILE_ENCODER
#include "lodepng.h"
#include "lodepng.cpp"

#include "bfbii.h"

#include "color_lib.h"

#define SIDE_Y 512
#define SIDE_X SIDE_Y


void yuv2rgb(double y, uint8_t cb, uint8_t cr, uint8_t *r, uint8_t *g, uint8_t *b)
{
	double r1, b1, g1;

	// Table with the 1's complement translation of chroma
	double c_table[32] = {
		-15, -14, -13, -12, -11, -10, -9, -8,
		-7,  -6,  -5,  -4,  -3,  -2,  -1, -0,
		 0,   1,   2,   3,   4,   5,   6,  7,
		 8,   9,   10,  11,  12,  13, 14,  15
	};

	// Y is 0 to 31
	y = (y / 31.0);

	// Populate color coefficients to NTSC REC 601
	struct color_coeff co = {0};
	pop_rec_const(REC_601, &co);

	// Normalize to MAX to -MAX for each signal. 14/15 roughly represents the standard
	// reduction constants (235/255) that would be applicable to these signals. Since the
	// table never goes below 1 or above 30 It's assumed this was intentional.
	// Note that V_MAX seems slightly too low for tg16, which is why it is slightly reduced.

	double pb = ((double) c_table[cb] / (c_table[31])) * co.u_max;
	double pr = ((double) c_table[cr] / (c_table[31])) * (co.v_max * 0.96);

	// Apply phase shift. This appears to vary a bit per system, but it definitely exists
	// intentionally as without it, the colors would be non-sensical. This was based on
	// a DUO-R belonging to Artemio.

	double phase_deg = -10.4; // in degrees
	double rad = deg_to_rad(phase_deg);
	double s = sin(rad);
	double c = cos(rad);
	double old_pb = pb;
	pb = pb * c - pr * s;

	rad = deg_to_rad(phase_deg);
	s = sin(rad);
	c = cos(rad);
	pr = pr * c + old_pb * s;

	// Reference NTSC YUV to RGB conversion.
	yuv_to_rgb(&co, y, pb, pr, &r1, &g1, &b1);

	// Desaturate the colors to the range the system purportedly used. This is very close on the
	// vectorscope.
	desaturate_rgb_linear(&co, 192.0/255.0, &r1, &g1, &b1);

	// Create 8 bit colors
	double m = 255.0;
	r1 *= m;
	g1 *= m;
	b1 *= m;

	// Clip if/when needed
	r1 = r1 > 255 ? 255 : r1 < 0 ? 0 : r1;
	g1 = g1 > 255 ? 255 : g1 < 0 ? 0 : g1;
	b1 = b1 > 255 ? 255 : b1 < 0 ? 0 : b1;

	*r = round(r1);
	*g = round(g1);
	*b = round(b1);
}


int main (int *argc, char *argv)
{
	uint8_t rgb_lut[512][3] = {0};

	uint8_t r, g, b;

	for (int32_t x = 0; x < 512; x++) {
		uint16_t index = ((index_lut[x][0] << 6) | (index_lut[x][1] << 3) | (index_lut[x][2])) & 0x1FF;
		int8_t y, u, v;
		yuv2rgb(pbpry_int_lut[x][2], pbpry_int_lut[x][0], pbpry_int_lut[x][1], &r, &g, &b);
		rgb_lut[index][0] = r;
		rgb_lut[index][1] = g;
		rgb_lut[index][2] = b;
	}

	unsigned char *raw_image = calloc(1024 * 1024 * 3, 1);

	for (int32_t x = 0; x < 512; x++) {
		printf("%3d: (%3d %3d %3d)\n", x,
			rgb_lut[x][0], rgb_lut[x][1], rgb_lut[x][2]
		);
		// Draw the palettes
		size_t ind = x * SIDE_X * 3;
		for (int y = 0; y < (SIDE_X * 3); y+=3) {
			raw_image[ind + y] = rgb_lut[x][0];
			raw_image[ind + y + 1] = rgb_lut[x][1];
			raw_image[ind + y + 2] = rgb_lut[x][2];
		}
	}

	// Make the files
	unsigned char *pngfile = NULL;
	size_t pngsize;

	lodepng_encode_memory(&pngfile, &pngsize, raw_image, SIDE_X, SIDE_Y, 2, 8);

	FILE *f = fopen("input1.png", "wb");
	fwrite(pngfile, 1, pngsize, f);
	fclose(f);

	free(raw_image);
	free(pngfile);

	// MiSTer Palette
	f = fopen("palette.tgp", "wb");
	for (int32_t x = 0; x < 512; x++) {
		uint16_t rgb[3] = {0,0,0};

		rgb[0] = rgb_lut[x][0];
		rgb[1] = rgb_lut[x][1];
		rgb[2] = rgb_lut[x][2];
		fwrite(rgb, 1, 3 * sizeof(uint16_t), f);
	}
	fclose(f);

	// Mednafen Palette
	f = fopen("palette.pal", "wb");
	for (int32_t x = 0; x < 512; x++) {
		uint8_t rgb[3] = {0,0,0};

		rgb[0] = rgb_lut[x][0];
		rgb[1] = rgb_lut[x][1];
		rgb[2] = rgb_lut[x][2];
		fwrite(rgb, 1, 3 * sizeof(uint8_t), f);
	}
	fclose(f);
}
