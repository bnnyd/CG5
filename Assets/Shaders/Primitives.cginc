void updateRayHit(inout RayHit bestHit, float3 pos, float t, float3 normal, Material material) {
	// update if the current intersection is at a valid distance and its distance is minimal 
	if ((t > 0) && ((t < bestHit.distance) || (isinf(bestHit.distance)))) {
		bestHit.position = pos;
		bestHit.distance = t;
		bestHit.normal = normal;
		bestHit.material = material;
	}
}

float4 planeRayIntersection(float3 orig, float3 dir, float3 planePoint, float3 planeNormal) {

	float cos_theta = dot(dir, planeNormal);

	if (cos_theta == 0) {
		// the plane lays on the same direction of the ray, the plane is not visible   
		return -1;
	}

	float t = -dot(orig - planePoint, planeNormal) / cos_theta;

	float3 position = orig + dir * t;
	return float4(position, t);
}

float solveTQuadraticEquation(float a, float b, float c){
	// check the discriminant delta and then solve the quadratic equation   
	float delta = b * b - 4 * a*c;
	float t;

	if (delta < 0) {
		// no intestections  
		return -1;

	}
		
	if (delta == 0) {
		t = -b / (2 * a);
	}
	else {
		t = min((-b - sqrt(delta)) / (2 * a), (-b + sqrt(delta)) / (2 * a));
		if (t <= 0) { // update the solution t if the smaller solution is negative     
			t = max((-b - sqrt(delta)) / (2 * a), (-b + sqrt(delta)) / (2 * a));
		}
	}
	return t;
}

// Checks for an intersection between a ray and a sphere
// The sphere center is given by sphere.xyz and its radius is sphere.w
void intersectSphere(Ray ray, inout RayHit bestHit, Material material, float4 sphere)
{
	float3 center = sphere.xyz;
	float radius = sphere.w;
	float3 dir = ray.direction;
	float3 orig = ray.origin;

	// solve the equation:
	// dir*dir t^2 + 2(orig - center)*dir t + (orig - center)(orig - center) - radius^2 = 0
	// call a,b,c the scalars that multiply t^2, t, 1 respectivelly

	float a = 1;
	float b = 2 * dot(orig - center, dir);
	float c = dot(orig - center, orig - center) - pow(radius, 2);

	float t = solveTQuadraticEquation(a, b, c);
	if (t < 0) {
		return;
    }
	float3 position = orig + dir * t;
	float3 normal = normalize(position - center);
	updateRayHit(bestHit, position, t, normal, material);
	return;
}

// Checks for an intersection between a ray and a plane
// The plane passes through point c and has a surface normal n 
void intersectPlane(Ray ray, inout RayHit bestHit, Material material, float3 c, float3 n)
{
	float3 dir = ray.direction;
	float3 orig = ray.origin;

	float4 intersection = planeRayIntersection(orig, dir, c, n);
	float3 position = intersection.xyz;
	float t = intersection.w;
	updateRayHit(bestHit, position, t, n, material);
	return;

}

Material calcCheckeredMaterial(Material m1, Material m2, float3 intersectionPoint, float3 planeNormal)
{
	float u;
	float v;

	if (planeNormal.x != 0){
		u = intersectionPoint.y;
		v = intersectionPoint.z;	
    } else if (planeNormal.y != 0){
		u = intersectionPoint.x;
		v = intersectionPoint.z;
	} else {
		u = intersectionPoint.x;
		v = intersectionPoint.y;
	}

	// xor: if u and v disagree is True:  
	if ( ((u - floor(u)) < 0.5) ^ ((v - floor(v)) < 0.5) )	
	{
		return m2;
	}
	else
	{
		return m1;
	}
}

// Checks for an intersection between a ray and a plane
// The plane passes through point c and has a surface normal n
// The material returned is either m1 or m2 in a way that creates a checkerboard pattern 
void intersectPlaneCheckered(Ray ray, inout RayHit bestHit, Material m1, Material m2, float3 c, float3 n)
{
	float3 dir = ray.direction;
	float3 orig = ray.origin;
	float3 planePoint = c;
	float3 planeNormal = normalize(n);

	float4 intersection = planeRayIntersection(orig, dir, planePoint, planeNormal);
	float3 p = intersection.xyz;
	float t = intersection.w;
	// Check with Binyamin if we check that ray hits plane on the normal side or the other
	Material material = calcCheckeredMaterial(m1, m2, p, planeNormal);
	updateRayHit(bestHit, p, t, planeNormal, material);
	return;
}


// Checks for an intersection between a ray and a triangle
// The triangle is defined by points a, b, c
void intersectTriangle(Ray ray, inout RayHit bestHit, Material material, float3 a, float3 b, float3 c)
{
	// find point p: meeting of the trinagle plane and the ray:  
	float3 dir = ray.direction;
	float3 orig = ray.origin;
	float3 planePoint = a;
	float3 planeNormal = normalize(cross(a - c, b - c));

	float4 intersection = planeRayIntersection(orig, dir, planePoint, planeNormal);
	float3 p = intersection.xyz;
	float t = intersection.w;
	if (t < 0) {
		// the ray doen't meet the plane in any valid point
		return;
	}
	// p is inside the triangle if all the the following products are non-negative:
	float prod1 = dot(cross(b - a, p - a), planeNormal);
	float prod2 = dot(cross(c - b, p - b), planeNormal);
	float prod3 = dot(cross(a - c, p - c), planeNormal);

	if ((prod1 < 0) || (prod2 < 0) || (prod3 < 0)) {
		return;
	}

	updateRayHit(bestHit, p, t, planeNormal, material);
	return;
}


// Checks for an intersection between a ray and a 2D circle
// The circle center is given by circle.xyz, its radius is circle.w and its orientation vector is n 
void intersectCircle(Ray ray, inout RayHit bestHit, Material material, float4 circle, float3 n)
{
	float3 dir = ray.direction;
	float3 orig = ray.origin;
	float3 planePoint = circle.xyz;
	float3 planeNormal = n;

	float4 intersection = planeRayIntersection(orig, dir, planePoint, planeNormal);
	float3 p = intersection.xyz;
	float t = intersection.w;
	if (t < 0) {
		// the ray doen't meet the plane in any valid point
		return;
	}
	if (distance(p, circle.xyz) > circle.w) {
		// the point is out of the circle
		return;
    }
	updateRayHit(bestHit, p, t, planeNormal, material);
	return;
}


// Checks for an intersection between a ray and a cylinder aligned with the Y axis
// The cylinder center is given by cylinder.xyz, its radius is cylinder.w and its height is h
void intersectCylinderY(Ray ray, inout RayHit bestHit, Material material, float4 cylinder, float h)
{
	float3 center = cylinder.xyz;
	float radius = cylinder.w;
	float height = h;
	float3 dir = ray.direction;
	float3 orig = ray.origin;

	float a = dir.x*dir.x + dir.z*dir.z;
	float b = 2 * ((orig.x - center.x)*dir.x + (orig.z - center.z)*dir.z);
	float c = pow(orig.x - center.x, 2) + pow(orig.z - center.z, 2) - pow(radius, 2);

	float t = solveTQuadraticEquation(a, b, c);
	if (t < 0) {
		return;
    }

	float3 position = orig + dir * t;

	if(position.y - center.y > height/2){
		float3 upCircleCenter = float3(center.x, center.y + height/2, center.z);
		float3 upNormal = float3(0,1,0);
		float4 upCircle = float4(upCircleCenter, radius);
		intersectCircle(ray, bestHit, material, upCircle, upNormal);

	} else if(center.y - position.y > height/2){
		float3 downCircleCenter = float3(center.x, center.y - height/2, center.z);
		float3 downNormal = float3(0,-1,0);
		float4 downCircle = float4(downCircleCenter, radius);
		intersectCircle(ray, bestHit, material, downCircle, downNormal);

	} else {
		float3 normal_temp = normalize(position - center);
		// project normal_temp on XZ plane:
		float3 normal = normalize(float3(normal_temp.x, 0, normal_temp.z));
		updateRayHit(bestHit, position, t, normal, material);
	}
	
	return;
}
