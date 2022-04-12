#pragma once

// Mini-library of color conversion functions for older broadcast standards
// Copyright Jamie Dickson, 2020.

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>


// Use POSIX M_PI if available
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

// RGB Coefficients per specification are fixed precision
// Fun trivia: CO_R + CO_G + CO_B == 1

#define CO_R_601 0.299
#define CO_B_601 0.114
#define CO_G_601 (1.0 - CO_R_601 - CO_B_601)

#define CO_R_709 0.2126
#define CO_B_709 0.0722
#define CO_G_709 (1.0 - CO_R_709 - CO_B_709)

#define CO_R_2020 0.2627
#define CO_B_2020 0.0593
#define CO_G_2020 (1.0 - CO_R_2020 - CO_B_2020)

typedef enum color_rec {
	REC_601 = 0,
	REC_709 = 1,
	REC_2020 = 2,
} COLOR_REC_T;

struct color_coeff {
	double r;
	double g;
	double b;
	double u_max;
	double v_max;
};

// Max and Min values as defined in BT 601 for Chroma. Luma may vary for NTSC-J.
#define MAX_IRE 130.8333
#define MIN_IRE 23.3025

// U V reduction factors for NTSC video
#define U_REDUCTION 0.492111
#define V_REDUCTION 0.877283

void pop_rec_const(COLOR_REC_T rec, struct color_coeff *co)
{
	memset(co, 0, sizeof(struct color_coeff));

	switch(rec) {
		case REC_709:
			co->r = CO_R_709;
			co->g = CO_G_709;
			co->b = CO_B_709;
		case REC_2020:
			co->r = CO_R_2020;
			co->g = CO_G_2020;
			co->b = CO_B_2020;
		default: // 601 as default
			co->r = CO_R_601;
			co->g = CO_G_601;
			co->b = CO_B_601;
	}

	co->u_max = 0.436;//(co->g + co->r);
	co->v_max = 0.615;//(co->g + co->b);
}

double ire_to_mv(double ire)
{
	return (ire * 7.14);
}

double mv_to_ire(double mv)
{
	return (mv / 7.14);
}

double deg_to_rad(double deg)
{
	return deg * (M_PI / 180.0);
}

double rad_to_deg(double rad)
{
	return rad * 180.0 / M_PI;
}

double uv_to_deg(double u, double v)
{
	double ang = atan2(v, u);

	return rad_to_deg(ang) + 180.0;
}

// Quantifies a normalized y value to ntsc-j standards
uint8_t y_quant_j(double y)
{
	int32_t y1 = round(y * 235.0);
	if (y1 < 0) y1 = 0;
	if (y1 > 255) y1 = 255;

	return y1;
}

// Quantifies a normalized y value to ntsc standards
uint8_t y_quant(double y)
{
	int32_t y1 = round((y * 219.0) + 16.0);
	if (y1 < 0) y1 = 0;
	if (y1 > 255) y1 = 255;

	return y1;
}

// Quantifies a normalized pb or pr value to rec 601 standards
uint8_t c_quant_601(double c)
{
	int32_t c1 = round((c * 219.0) + 16.0);
	if (c1 < 0) c1 = 0;
	if (c1 > 255) c1 = 255;

	return c1;
}

// Quantifies a normalized pb or pr value to rec 709 standards
uint8_t c_quant_709(double c)
{
	int32_t c1 = round((c * 224.0) + 16.0);
	if (c1 < 0) c1 = 0;
	if (c1 > 255) c1 = 255;

	return c1;
}

// Standard split chroma to reduced u and v
static void c_to_uv(double c, double *u, double *v)
{
	*v = sin(c);
	*u = cos(c);
}

// Standard RGB to luma generation for most color spaces
static double rgb_to_y(struct color_coeff *co, double r, double g, double b)
{
	return (co->r * r) + (co->g * g) + (co->b * b);
}

// Convert standard RGB to un-reduced YUV
static void rgb_to_yuv(struct color_coeff *co, double r, double g, double b, double *y, double *u, double *v)
{

	*y = rgb_to_y(co, r, g, b);

	*u = (co->u_max * b) - (co->r * r) - (co->g * g);
	*v = (co->v_max * r) - (co->g * g) - (co->b * b);
}

// Convert un-reduced YUV to RGB
static void yuv_to_rgb(struct color_coeff *co, double y, double u, double v, double *r, double *g, double *b)
{
	*r = y + v * ((1.0 - co->r) / co->v_max);
	*g = y - u * ((co->b * (1.0 - co->b)) / (co->u_max * co->g)) - v * ((co->r * (1.0 - co->r)) / (co->v_max * co->g));
	*b = y + u * ((1.0 - co->b) / co->u_max);
}

void rgb_to_ycbcr(struct color_coeff *co, double r, double g, double b, double *y, double *cb, double *cr)
{

	*y = rgb_to_y(co, r, g, b);

	*cb = b - (co->r * r) - (co->g * g);
	*cr = r - (co->g * g) - (co->b * b);
}

void ycbcr_to_rgb(struct color_coeff *co, double y, double cb, double cr, double *r, double *g, double *b)
{
	*r = y + cr * ((1.0 - co->r));
	*g = y - cb * ((co->b * (1.0 - co->b)) / (co->g)) - cr * ((co->r * (1.0 - co->r)) / (co->g));
	*b = y + cb * ((1.0 - co->b));
}

void desaturate_rgb(double factor, double *r, double *g, double *b)
{
	double L = 0.3 * *r + 0.6 * *g + 0.1 * *b;
	*r = *r + factor * (L - *r);
	*g = *g + factor * (L - *g);
	*b = *b + factor * (L - *b);
}
void desaturate_rgb_linear(struct color_coeff *co, double factor, double *r, double *g, double *b)
{
	double y = rgb_to_y(co, *r, *g, *b);
	*r = (*r * factor) + (y * (1.0 - factor));
	*g = (*g * factor) + (y * (1.0 - factor));
	*b = (*b * factor) + (y * (1.0 - factor));
}
// // Convert RGB to YCbCr
// void rgb_to_ycbcr(struct color_coeff *co, double r, double g, double b, double *y, double *cb, double *cr)
// {
// 	*y = rgb_to_y(co, r, g, b);
// 	*cb = 0.5 * ((b - *y) / (1.0 - co->b));
// 	*cr = 0.5 * ((r - *y) / (1.0 - co->r));
// }

// // This does not assume a bit size, so scaling and reductions should be done afterwards
// void ycbcr_to_rgb(struct color_coeff *co, double y, double cb, double cr, double *r, double *g, double *b)
// {
// 	*r = y + (co->v_max * 2.0) * cr;
// 	*g = y - ((co->u_max * 2.0) * co->b) / (co->g * cb) - ((co->v_max * 2.0) * co->r) / (co->g * cr);
// 	*b = y + (co->u_max * 2.0) * cb;
// }

// Creates a standard rec 601 U and V reduced value for composite video from (B-Y) and (R-Y)
static void uv_reduce(double *u, double *v)
{
	// For compatible with 50's era television recievers, it was anecodtally discovered that the
	// color difference signals must not have a possible excursion of greater than 33.3% above white
	// or below black. This means a maximum positive excursion of 100IRE + ((1/3) * 92.5 IRE) =
	// 130.8025 IRE and a minimum of 7.5IRE - ((1/3) * 92.5) = 23.3025IRE. The goal is to make the
	// yellow and cyan bars of 75% color bars equal to 100 IRE.
	*u = U_REDUCTION * *u;
	*v = V_REDUCTION * *v;
}

static void uv_to_iq(double u, double v, double *i, double *q)
{
	// IQ is effectively UV rotated about 33 degrees for ancient bandwidth and fidelity reasons,
	// with some attempts made to align the signal to have more fidelity with skin tones. Contrary
	// to popular belief, IQ is not used for most modern-ish composite signals, but rather UV is
	// instead.
	*i = (u * sin(deg_to_rad(33.0)) + (v * cos(deg_to_rad(33.0))));
	*q = (u * cos(deg_to_rad(33.0)) + (v * sin(deg_to_rad(33.0))));
}

// Angle + Vector to X Y coordinates
void av_to_xy(double angle, double vector, double *x, double *y)
{
	*x = vector * cos(deg_to_rad(angle));
	*y = vector * sin(deg_to_rad(angle));
}

// Generate the maximum possible vector for a given angle
double angle_max_vector(struct color_coeff *co, double angle)
{
	double u, v;

	u = (255.0 * co->u_max) * cos(deg_to_rad(angle));
	v = (255.0 * co->v_max) * sin(deg_to_rad(angle));

	u *= U_REDUCTION;
	v *= V_REDUCTION;

	return sqrt((u * u) + (v * v));
}

// Returns the reduced vector after applying NTSC U and V reduction for a given angle.
double av_reduce(double angle, double vector)
{
	double v1 = 0;
	double x, y;
	av_to_xy(angle, vector, &x, &y);

	uv_reduce(&x, &y);

	v1 = sqrt((x * x) + (y * y));

	return v1;
}