/*
 * iapws - IAPWS formulations for the properties of water and steam
 * Copyright (C) 2022 Jonathan Debove
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef IAPWS_COEF_H
#define IAPWS_COEF_H

#include "pow.h"

struct Ni {
	int i;
	double N;
};

struct Nij {
	int i;
	int j;
	double N;
};

#define COEF_ITERATE1(coef, len, x, loop_body) do {	\
	int _ind_ = 0;					\
	double _xi_ = POWINT((x), (coef)[_ind_].i);	\
	double _ans_ = (coef)[_ind_].N * _xi_;		\
	loop_body;					\
	for (_ind_ = 1; _ind_ < (len); ++_ind_) {	\
		_xi_ = (coef)[_ind_].i >= (coef)[_ind_ - 1].i ? _xi_ * POWINT((x), (coef)[_ind_].i - (coef)[_ind_ - 1].i) : POWINT((x), (coef)[_ind_].i);	\
		_ans_ = (coef)[_ind_].N * _xi_;		\
		loop_body;				\
	}						\
} while(0)

#define COEF_ITERATE2(coef, len, x, y, loop_body) do {	\
	int _ind_ = 0;					\
	double _xi_ = POWINT((x), (coef)[_ind_].i);	\
	double _yj_ = POWINT((y), (coef)[_ind_].j);	\
	double _ans_ = (coef)[_ind_].N * _xi_ * _yj_;	\
	loop_body;					\
	for (_ind_ = 1; _ind_ < (len); ++_ind_) {	\
		_xi_ = (coef)[_ind_].i >= (coef)[_ind_ - 1].i ? _xi_ * POWINT((x), (coef)[_ind_].i - (coef)[_ind_ - 1].i) : POWINT((x), (coef)[_ind_].i);	\
		_yj_ = (coef)[_ind_].j >= (coef)[_ind_ - 1].j ? _yj_ * POWINT((y), (coef)[_ind_].j - (coef)[_ind_ - 1].j) : POWINT((y), (coef)[_ind_].j);	\
		_ans_ = (coef)[_ind_].N * _xi_ * _yj_;	\
		loop_body;				\
	}						\
} while(0)

#endif
