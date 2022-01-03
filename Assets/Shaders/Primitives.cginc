void updateRayHit(inout RayHit bestHit, float3 pos, float t, float3 normal, Material material){
    // update if the current intersection is at a valid distance and its distance is minimal 
    if ((t > 0) && (t < bestHit.distance)){
        bestHit.position = pos;
        bestHit.distance = t;
        bestHit.normal = normal;
        bestHit.material = material;
        }
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

        float t = (- b - sqrt(delta))/(2*a);

        // update the solution t if the smaller solution is negative
        if ((delta > 0) && (t <= 0)){
            t = (- b + sqrt(delta))/(2*a);
        }

        float3 position = orig + dir*t;
        float3 normal = normalize(bestHit.position - center);
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

    float cos_theta = dot(dir, n);
    if(cos_theta == 0){
        // the plane lays on the same direction of the ray, the plane is not visible 
        return;
    }
    float t = -dot(orig - c, n)/cos_theta;
    float3 position = orig + dir*t;

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
    // Your implementation
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
