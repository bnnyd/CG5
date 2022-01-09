void updateRayHit(inout RayHit bestHit, float3 pos, float t, float3 normal, Material material){
    // update if the current intersection is at a valid distance and its distance is minimal 
    if ((t > 0) && ((t < bestHit.distance) || (isinf(bestHit.distance)))){
        bestHit.position = pos;
        bestHit.distance = t;
        bestHit.normal = normal;
        bestHit.material = material;
        }
}

float4 planeRayIntersection(float3 orig, float3 dir, float3 planePoint, float3 planeNormal){
    
    float cos_theta = dot(dir, planeNormal);

    if(cos_theta == 0){
        // the plane lays on the same direction of the ray, the plane is not visible   
        return -1;
    }

    float t = -dot(orig - planePoint, planeNormal)/cos_theta;
   
    float3 position = orig + dir*t;
    return float4(position, t);
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
    float b = 2*dot(orig - center, dir);
    float c = dot(orig - center, orig - center) - pow(radius, 2);
    
    // check the discriminant delta and then solve the quadratic equation

    float delta = b*b - 4*a*c;

    if(delta < 0){
        // no intestections
        return;

    } else {

        float t;
        if(delta == 0) {
            t = -b/(2*a);
        } else{
            t = min((- b - sqrt(delta))/(2*a), (- b + sqrt(delta))/(2*a));
            if (t <= 0){ // update the solution t if the smaller solution is negative  
                t = max((- b - sqrt(delta))/(2*a), (- b + sqrt(delta))/(2*a));
            }
        }
        
        float3 position = orig + dir*t;
        float3 normal = normalize(position - center);
        updateRayHit(bestHit, position, t, normal, material);
        return;      
    }
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

// Checks for an intersection between a ray and a plane
// The plane passes through point c and has a surface normal n
// The material returned is either m1 or m2 in a way that creates a checkerboard pattern 
void intersectPlaneCheckered(Ray ray, inout RayHit bestHit, Material m1, Material m2, float3 c, float3 n)
{
    // Your implementation
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
    if(t < 0){
        // the ray doen't meet the plane in any valid point
        return;
    }
    // p is inside the triangle if all the the following products are non-negative:
    float prod1 = dot(cross(b-a,p-a),planeNormal);
    float prod2 = dot(cross(c-b,p-b),planeNormal);
    float prod3 = dot(cross(a-c,p-c),planeNormal);

    if((prod1 < 0) || (prod2 < 0) || (prod3 < 0)){
        return;
    }

    updateRayHit(bestHit, p, t, planeNormal, material);
    return;
}


// Checks for an intersection between a ray and a 2D circle
// The circle center is given by circle.xyz, its radius is circle.w and its orientation vector is n 
void intersectCircle(Ray ray, inout RayHit bestHit, Material material, float4 circle, float3 n)
{
    // Your implementation
}


// Checks for an intersection between a ray and a cylinder aligned with the Y axis
// The cylinder center is given by cylinder.xyz, its radius is cylinder.w and its height is h
void intersectCylinderY(Ray ray, inout RayHit bestHit, Material material, float4 cylinder, float h)
{
    // Your implementation
}
