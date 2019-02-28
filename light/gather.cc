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

#include <light/light.hh>
#include <light/entities.hh>
#include <light/gather.hh>
#include <light/ltface.hh>

#include <common/mathlib.hh>

#include <cassert>

using namespace std;


/* Dirtmapping borrowed from q3map2, originally by RaP7oR */
rayGenFunc dirtRayGen = 0;
rayGenFunc domeRayGen = 0;
std::vector<qvec3f> dirtVectors;

static inline float crandom() { return 1 - Random() * 2; }
static inline float lerp(float a, float b, float mix) { return a + (b - a)*mix; }

// from q3map2
void
GetUpRtVecs(const vec3_t normal, vec3_t myUp, vec3_t myRt)
{
	/* check if the normal is aligned to the world-up */
	if (normal[0] == 0.0f && normal[1] == 0.0f) {
		if (normal[2] == 1.0f) {
			VectorSet(myRt, 1.0f, 0.0f, 0.0f);
			VectorSet(myUp, 0.0f, 1.0f, 0.0f);
		}
		else if (normal[2] == -1.0f) {
			VectorSet(myRt, -1.0f, 0.0f, 0.0f);
			VectorSet(myUp, 0.0f, 1.0f, 0.0f);
		}
	}
	else {
		vec3_t worldUp;
		VectorSet(worldUp, 0.0f, 0.0f, 1.0f);
		CrossProduct(normal, worldUp, myRt);
		VectorNormalize(myRt);
		CrossProduct(myRt, normal, myUp);
		VectorNormalize(myUp);
	}
}

// from q3map2
void
TransformToTangentSpace(const vec3_t normal, const vec3_t myUp, const vec3_t myRt, const qvec3f inputvec, vec3_t outputvec)
{
	for (int i = 0; i < 3; i++)
		outputvec[i] = myRt[i] * inputvec[0] + myUp[i] * inputvec[1] + normal[i] * inputvec[2];
}
void
TransformToTangentSpace(const vec3_t normal, const vec3_t up, const vec3_t rt, const qvec3f &inputvec, qvec3f &outputvec)
{
	qvec3f temp;
	for (int i = 0; i < 3; i++) {
		temp[i] = rt[i] * inputvec[0] + up[i] * inputvec[1] + normal[i] * inputvec[2];
	}
	outputvec = temp;
}



/* ======================================================================== 
 *
 *		RAY GENERATORS
 *
 * ======================================================================== */

/*
 * ============
 * RandomRay
 * 
 * generate a random unit vector within 'maxpitch' degrees of (0 0 1)
 * z is biased by cos(declination) to avoid polar clustering around z+
 * ============
 */
static qvec3f
RandomRay(const float maxpitch)
{
	float ang, z, k, cmp;

	ang = 2 * Q_PI * Random();	// yaw
	cmp = cos(maxpitch * Q_PI / 180.0);	// declination
	z = lerp(cmp, 1.0, Random());	// eliminate bunching up at the poles
	k = sqrt(1 - z * z);
	return qvec3f(k * cos(ang), k * sin(ang), z);
}


/*
========================================================================
Random Ray Generation
- use simple random sampling
- rays are technically uniformly distributed but can appear clustered, which
	manifests visually as dirt shadows being cast more strongly in certain
	directions
*/

/*
 * ============
 * GenRaysRandom
 * ============
 */
static int
GenRaysRandom(const int samples, const float maxpitch, std::vector<qvec3f> &rays)
{
	rays.clear();
	rays.reserve(samples);

	for (int i = 0; i < samples; i++) {
		rays.push_back(RandomRay(maxpitch));
	}

	return samples;
}





/*
========================================================================
Elevation/Angle Ray Generation
- the original, default dirt ray generator taken from q3map2
- does not compensate properly for rays bunching up near the pole, which
	manifests visually as dirt gamma being a little higher than it is
*/

static int
GenRaysAngle(const int desiredSamples, const float maxpitch, std::vector<qvec3f> &rays)
{
	int numRays = 0;

	/* lunaran - pick angle and elevation steps that yield close to 
	 * desiredSamples with roughly even spacing in latitude and longitude */
	float r = sqrt(desiredSamples * 360.0f / maxpitch);
	const int pitchSteps = qmax(2, ceil(desiredSamples / r));
	const int angleSteps = desiredSamples / pitchSteps;
	const float angleStep = (float)DEG2RAD(360.0f / angleSteps);
	const float elevationStep = (float)DEG2RAD(maxpitch / pitchSteps);
	//logprint("%i samp, %f pitch: %i x %i = %i\n", desiredSamples, maxpitch, angleSteps, pitchSteps, angleSteps * pitchSteps);

	rays.clear();
	rays.reserve(pitchSteps * angleSteps);

	/* iterate angle */
	float angle = 0.0f;
	for (int i = 0; i < angleSteps; i++, angle += angleStep) {
		/* iterate elevation */
		float elevation = elevationStep * 0.5f;
		for (int j = 0; j < pitchSteps; j++, elevation += elevationStep) {
			rays.emplace_back(sin(elevation) * cos(angle), sin(elevation) * sin(angle), cos(elevation));
			numRays++;
		}
	}
	return numRays;
}




/*
========================================================================
Poisson Disc Ray Generation
- rays are randomly generated, but rejected if they are not within a distance
	threshold of their neighbors
- ray are not uniformly distributed, but appear at least evenly distributed,
	and appearances are what matter here
- significantly slower to generate large amounts of rays than halton
- usually can't generate the exact number of rays requested but it sure tries
- https://www.jasondavies.com/poisson-disc/
*/

static float
distsqr(const vec3_t a, const vec3_t b)
{
	vec3_t c;
	VectorSubtract(a, b, c);
	return c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
}

static qboolean
CheckRayAnnulus(const vec3_t c, const vec3_t r, const float dist)
{
	float dsq = distsqr(c, r);
	float inner = dist * dist;
	float outer = 4 * inner;
	//logprint("  distance: %f\n", sqrtf(dsq));
	if (dsq < inner) {
		return false;
	}
	return (dsq < outer);
}

/*
 * ============
 * GenRaysPoisson
 * ============
 */
static int
GenRaysPoisson(const int desiredSamples, const float maxpitch, vec3_t** rays)
{
	int activeRay, activeRayIdx, numRays, numActiveRays;
	int i, j;
	float cmp;

	cmp = cos(DEG2RAD(maxpitch));	/* max declination dot */

	/* guess at a good average distance between rays in order to fit the desired number
	 of rays into a cone maxpitch*2 degrees wide */
	float angSize = maxpitch / (40 * sqrt(desiredSamples) * (1 - 1 / desiredSamples));

	*rays = (vec3_t*)calloc((desiredSamples + 16), sizeof(vec3_t));
	int* active = (int*)calloc((desiredSamples + 16), sizeof(int));

	/* starting ray */
	VectorSet((*rays)[0], 0, 0, 1);
	numRays = 1;
	active[0] = 0;
	numActiveRays = 1;

	while (numActiveRays > 0) {
		activeRayIdx = (int)floor(Random() * numActiveRays);
		activeRay = active[activeRayIdx];

		if (activeRay == -1) continue;

		vec3_t up, rt;
		GetUpRtVecs((*rays)[activeRay], up, rt);
		qboolean found = false;
		for (i = 0; i < 10; i++) {
			vec3_t candidate;
			/* pick random rays until one is within the annulus */
			for (j = 0; j < 32; j++) {
				VectorCopy((*rays)[activeRay], candidate);
				VectorMA(candidate, angSize * crandom(), up, candidate);
				VectorMA(candidate, angSize * crandom(), rt, candidate);
				VectorNormalize(candidate);

				if (candidate[2] < cmp) {
					continue;	/* outside maxpitch cone */
				}
				if (CheckRayAnnulus(candidate, (*rays)[activeRay], angSize)) {
					break;		/* found a suitable new ray  */
				}
			}
			if (j == 32) {
				//logprint("%s: can't find candidate for ray %i (%1.5f %1.5f %1.5f)", __func__, activeRay, rays[activeRay][0], rays[activeRay][1], rays[activeRay][2]);
				continue;
			}

			/* test it against all rays */
			qboolean passed = true;
			for (j = 0; j < numRays; j++) {
				if (distsqr(candidate, (*rays)[j]) < angSize*angSize) {
					passed = false;
					break;
				}
			}
			/* if it clears all checks, add it to the active list */
			if (passed) {
				found = true;
				VectorCopy(candidate, (*rays)[numRays]);
				active[numActiveRays] = numRays;
				numRays++;
				numActiveRays++;
				break;
			}
		}
		if (!found) {
			/* all failed, the space around the starting ray is deemed 'full'
			 * remove that ray from the active list */
			for (i = activeRayIdx; i < numActiveRays; i++) {
				active[i] = active[i + 1];
			}
			numActiveRays--;
		}
		if (numRays == (desiredSamples + 15)) {
			//logprint("%s: sample overflow (%i)\n", __func__, numRays);
			free(active);
			return numRays;
		}
	}
	//logprint("Targeted %i samples, sphere-packing yielded %i\n", desiredSamples, numRays);
	free(active);
	return numRays;
}

/*
static void
PoissonTest_Ang(const int desiredSamples, const float maxpitch)
{
	float guess = maxpitch / (39 * sqrt(desiredSamples) * (1 - 1 / desiredSamples));
	float angSize = guess;
	float angStep = angSize * 0.5f;

	vec3_t* raysTemp;

	for (int i = 0; i < 10; i++)
	{
		int actualSamples = GenRaysPoisson(desiredSamples, maxpitch, &raysTemp);
		if (actualSamples == desiredSamples)
			break;

		if (actualSamples < desiredSamples)
			angSize -= angStep;
		else
			angSize += angStep;

		angStep /= 2;
	}
	logprint("%i samp, %1.1f cone, %f guess, %f ang\n", desiredSamples, maxpitch, guess, angSize);
}
*/

void
GatherTest()
{
	std::vector<qvec3f> raysTemp;
	for (int pitch = 90; pitch > 0; pitch -= 10)
	{
		for (int samp = 12; samp <= 384; samp *= 2)
		{
			GenRaysAngle(samp, pitch, raysTemp);
		}
	}
}

/*
========================================================================
Halton Sequence Ray Generation
- generate pseudorandom uniform rays that don't	appear to cluster
- https://www.slideshare.net/luk036/n-sphere

- use a non-zero starting index for a 'randomized'-per-sample application
*/

static float
Halton(int idx, int base)
{
	int i;
	float f, out;

	f = 1;
	out = 0;
	i = idx;
	while (i > 0) {
		f /= base;
		out += f * (i % base);
		i = i / base;
	}
	return out;
}

/*
 * ============
 * GenRaysHalton
 * ============
*/
static int
GenRaysHalton(const int samples, const float maxpitch, std::vector<qvec3f> &rays)
{
	int i;
	float ang, z, k, cmp;

	rays.clear();
	rays.reserve(samples);

	for (i = 0; i < samples; i++) {
		ang = 2 * Q_PI * Halton(i, 2);
		cmp = cos(maxpitch * Q_PI / 180.0);
		z = lerp(cmp, 1.0, Halton(i, 3));
		k = sqrt(1 - z * z);
		rays.emplace_back(k * cos(ang), k * sin(ang), z);
	}

	return samples;
}





int
GenRaysCone(rayGenFunc gen, const int samples, const vec3_t center, const float conesize, std::vector<qvec3f> &rays)
{
	gen(samples, conesize, rays);

	// align rays to center
	vec3_t up, right, fwd;
	VectorCopy(center, fwd);
	VectorScale(fwd, -1, fwd);
	GetUpRtVecs(fwd, up, right);

	for (int i = 0; i < samples; i++) {
		TransformToTangentSpace(fwd, up, right, rays[i], rays[i]);
	}

	return rays.size();
}

int
GenRays(rayGenFunc gen, const int samples, std::vector<qvec3f> &rays)
{
	gen(samples, 180, rays);
	return rays.size();
}



/* ======================================================================== */



/* ======================================================================== */

static inline float fraction(float min, float val, float max) {
	if (val >= max) return 1.0;
	if (val <= min) return 0.0;

	return (val - min) / (max - min);
}


/*
 * ============
 * Dirt_GetScaleFactor
 *
 * returns scale factor for dirt/ambient occlusion
 * ============
 */
vec_t
Dirt_GetScaleFactor(const globalconfig_t &cfg, vec_t occlusion, const light_t *entity, const vec_t entitydist, const lightsurf_t *surf)
{
	vec_t light_dirtgain = cfg.dirtGain.floatValue();
	vec_t light_dirtscale = cfg.dirtScale.floatValue();
	bool usedirt;

	/* is dirt processing disabled entirely? */
	if (!dirt_in_use)
		return 1.0f;
	if (surf && surf->nodirt)
		return 1.0f;

	/* should this light be affected by dirt? */
	if (entity) {
		if (entity->dirt.intValue() == -1) {
			usedirt = false;
		}
		else if (entity->dirt.intValue() == 1) {
			usedirt = true;
		}
		else {
			usedirt = cfg.globalDirt.boolValue();
		}
	}
	else {
		/* no entity is provided, assume the caller wants dirt */
		usedirt = true;
	}

	/* if not, quit */
	if (!usedirt)
		return 1.0;

	/* override the global scale and gain values with the light-specific
	   values, if present */
	if (entity) {
		if (entity->dirtgain.floatValue())
			light_dirtgain = entity->dirtgain.floatValue();
		if (entity->dirtscale.floatValue())
			light_dirtscale = entity->dirtscale.floatValue();
	}

	/* early out */
	if (occlusion <= 0.0f) {
		return 1.0f;
	}

	/* apply gain */
	float outDirt = pow(occlusion, light_dirtgain);
	if (outDirt > 1.0f) {
		outDirt = 1.0f;
	}

	/* apply scale */
	outDirt *= light_dirtscale;
	if (outDirt > 1.0f) {
		outDirt = 1.0f;
	}

	/* lerp based on distance to light */
	if (entity) {
		// From 0 to _dirt_off_radius units, no dirt.
		// From _dirt_off_radius to _dirt_on_radius, the dirt linearly ramps from 0 to full, and after _dirt_on_radius, it's full dirt.

		if (entity->dirt_on_radius.isChanged()
			&& entity->dirt_off_radius.isChanged()) {

			const float onRadius = entity->dirt_on_radius.floatValue();
			const float offRadius = entity->dirt_off_radius.floatValue();

			if (entitydist < offRadius) {
				outDirt = 0.0;
			}
			else if (entitydist >= offRadius && entitydist < onRadius) {
				const float frac = fraction(offRadius, entitydist, onRadius);
				outDirt = frac * outDirt;
			}
		}
	}

	/* return to sender */
	return 1.0f - outDirt;
}

/*
 * =============
 * Gather_SetupDomeGen
 * =============
 */
void
Gather_SetupDomeGen(const globalconfig_t &cfg)
{
	if (!domeRayGen) {
		switch (cfg.dirtMode.intValue()) {
		case 2:
			domeRayGen = GenRaysAngle;
			logprint("Using elevation/angle dome ray generation\n");
			break;
		case 1:
			domeRayGen = GenRaysRandom;
			logprint("Using random dome ray generation\n");
			break;
		case 0:
		default:
			domeRayGen = GenRaysHalton;
			logprint("Using Halton dome ray generation\n");
		}
	}
}

/*
 * ============
 * SetupDirt
 *
 * sets up dirtmap (ambient occlusion)
 * ============
 */
void SetupDirt(globalconfig_t &cfg)
{
	// check if needed
	if (!cfg.globalDirt.boolValue()
		&& cfg.globalDirt.isLocked()) {
		// HACK: "-dirt 0" disables all dirtmapping even if we would otherwise use it.
		dirt_in_use = false;
		return;
	}

	if (cfg.globalDirt.boolValue()
		|| cfg.minlightDirt.boolValue()
		|| cfg.sunlight_dirt.boolValue()
		|| cfg.sunlight2_dirt.boolValue()) {
		dirt_in_use = true;
	}

	if (!dirt_in_use) {
		// check entities, maybe only a few lights use it
		for (const auto &light : GetLights()) {
			if (light.dirt.boolValue()) {
				dirt_in_use = true;
				break;
			}
		}
	}

	if (!dirt_in_use) {
		// dirtmapping is not used by this map. 
		return;
	}

	/* note it */
	logprint("--- SetupDirt ---\n");

	/* clamp dirtAngle */
	if (cfg.dirtAngle.floatValue() <= 1.0f) {
		cfg.dirtAngle.setFloatValueLocked(1.0f); // FIXME: add clamping API
	}
	if (cfg.dirtAngle.floatValue() >= 90.0f) {
		cfg.dirtAngle.setFloatValueLocked(90.0f);
	}

	/* 
	lunaran - user specified ray generator
	Halton is the new default, because it yields the best-looking results and does it the fastest
	*/
	switch (cfg.dirtMode.intValue()) {
		/*
	case 3:
		dirtRayGen = GenRaysPoisson;
		logprint("Using Poisson dirt ray generation\n");
		break;
		*/
	case 2:
		dirtRayGen = GenRaysAngle;
		logprint("Using elevation/angle dirt ray generation\n");
		break;
	case 1:
		dirtRayGen = GenRaysRandom;
		logprint("Using random dirt ray generation\n");
		break;
	case 0:
	default:
		dirtRayGen = GenRaysHalton;
		logprint("Using Halton dirt ray generation\n");
	}

	dirtRayGen(cfg.dirtSamples.intValue(), cfg.dirtAngle.floatValue(), dirtVectors);

	/* emit some statistics */
	logprint("%9d dirtmap vectors\n", dirtVectors.size());
}

static qvec3f
GetDirtVector(const globalconfig_t &cfg, int i)
{
	Q_assert(i < dirtVectors.size());

	if (cfg.dirtMode.intValue() == 1) {
		return RandomRay(cfg.dirtAngle.floatValue());
	}
	else {
		return dirtVectors[i];
	}
}

float
DirtAtPoint(const globalconfig_t &cfg, raystream_t *rs, const vec3_t point, const vec3_t normal, const modelinfo_t *selfshadow)
{
	if (!dirt_in_use) {
		return 0.0f;
	}

	vec3_t myUp, myRt;
	float occlusion = 0;

	// this stuff is just per-point

	GetUpRtVecs(normal, myUp, myRt);

	rs->clearPushedRays();

	for (int j = 0; j < dirtVectors.size(); j++) {
		// fill in input buffers
		qvec3f dirtvec = GetDirtVector(cfg, j);
		vec3_t dir;
		TransformToTangentSpace(normal, myUp, myRt, dirtvec, dir);
		rs->pushRay(j, point, dir, cfg.dirtDepth.floatValue(), selfshadow);
	}

	Q_assert(rs->numPushedRays() == dirtVectors.size());

	// trace the batch
	rs->tracePushedRaysIntersection();

	// accumulate hitdists
	for (int j = 0; j < dirtVectors.size(); j++) {
		if (rs->getPushedRayHitType(j) == hittype_t::SOLID) {
			const float dist = rs->getPushedRayHitDist(j);
			occlusion += qmin(cfg.dirtDepth.floatValue(), dist);
		}
		else {
			occlusion += cfg.dirtDepth.floatValue();
		}
	}

	// process the results.
	const vec_t avgHitdist = occlusion / dirtVectors.size();
	occlusion = 1 - (avgHitdist / cfg.dirtDepth.floatValue());
	return occlusion;
}






/*
 * ============
 * LightFace_CalculateDirt
 * ============
 */
void
LightFace_CalculateDirt(lightsurf_t *lightsurf)
{
	const globalconfig_t &cfg = *lightsurf->cfg;

	Q_assert(dirt_in_use);

	// batch implementation:

	vec3_t *myUps = (vec3_t *)calloc(lightsurf->numpoints, sizeof(vec3_t));
	vec3_t *myRts = (vec3_t *)calloc(lightsurf->numpoints, sizeof(vec3_t));

	// init
	for (int i = 0; i < lightsurf->numpoints; i++) {
		lightsurf->occlusion[i] = 0;
	}

	// this stuff is just per-point
	for (int i = 0; i < lightsurf->numpoints; i++) {
		GetUpRtVecs(lightsurf->normals[i], myUps[i], myRts[i]);
	}

	for (int j = 0; j < dirtVectors.size(); j++) {
		raystream_t *rs = lightsurf->stream;
		rs->clearPushedRays();

		// fill in input buffers

		for (int i = 0; i < lightsurf->numpoints; i++) {
			if (lightsurf->occluded[i])
				continue;

			qvec3f dirtvec = GetDirtVector(cfg, j);
			vec3_t dir;
			TransformToTangentSpace(lightsurf->normals[i], myUps[i], myRts[i], dirtvec, dir);
			rs->pushRay(i, lightsurf->points[i], dir, cfg.dirtDepth.floatValue(), lightsurf->modelinfo);
		}

		// trace the batch
		rs->tracePushedRaysIntersection();

		// accumulate hitdists
		for (int k = 0; k < rs->numPushedRays(); k++) {
			const int i = rs->getPushedRayPointIndex(k);
			if (rs->getPushedRayHitType(k) == hittype_t::SOLID) {
				float dist = rs->getPushedRayHitDist(k);
				lightsurf->occlusion[i] += qmin(cfg.dirtDepth.floatValue(), dist);
			}
			else {
				lightsurf->occlusion[i] += cfg.dirtDepth.floatValue();
			}
		}
	}

	// process the results.
	for (int i = 0; i < lightsurf->numpoints; i++) {
		vec_t avgHitdist = lightsurf->occlusion[i] / (float)dirtVectors.size();
		lightsurf->occlusion[i] = 1 - (avgHitdist / cfg.dirtDepth.floatValue());
	}

	free(myUps);
	free(myRts);
}

