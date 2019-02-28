/*  Copyright (C) 1996-1997  Id Software, Inc.

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

	See file, 'COPYING', for details.
*/

#ifndef __LIGHT_GATHER_H__
#define __LIGHT_GATHER_H__

const float MAX_SKY_DIST = 65536.0f;

typedef int(*rayGenFunc)(const int samples, const float maxpitch, std::vector<qvec3f> &rays);

extern rayGenFunc dirtRayGen;
extern rayGenFunc domeRayGen;

void SetupDirt(globalconfig_t &cfg);
int GenRaysCone(rayGenFunc gen, const int samples, const vec3_t center, const float conesize, std::vector<qvec3f>& rays);
int GenRays(rayGenFunc gen, const int samples, std::vector<qvec3f>& rays);
vec_t Dirt_GetScaleFactor(const globalconfig_t &cfg, vec_t occlusion, const light_t *entity, const vec_t entitydist, const lightsurf_t *surf);
void Gather_SetupDomeGen(const globalconfig_t & cfg);
void GetUpRtVecs(const vec3_t normal, vec3_t myUp, vec3_t myRt);
void TransformToTangentSpace(const vec3_t normal, const vec3_t myUp, const vec3_t myRt, const qvec3f inputvec, vec3_t outputvec);
void TransformToTangentSpace(const vec3_t normal, const vec3_t up, const vec3_t rt, const qvec3f & inputvec, qvec3f & outputvec);
void LightFace_CalculateDirt(lightsurf_t *lightsurf);
void GatherTest();


#endif /* __LIGHT_GATHER_H__ */

