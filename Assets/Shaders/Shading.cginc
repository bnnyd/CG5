//static const float SPECULARITIY = 0.4f;
// Implements an adjusted version of the Blinn-Phong lighting model 
float3 blinnPhong(float3 n, float3 v, float3 l, float shininess, float3 albedo)
{
    float3 SPECULARITIY = 0.4f;
	float3 diffuse = max(0, dot(n, l)) * albedo;
	float3 h = normalize(l + v);
	float3 specular = pow(max(0, dot(n, h)), shininess) * SPECULARITIY;
	return diffuse + specular;
}

// Reflects the given ray from the given hit point
void reflectRay(inout Ray ray, RayHit hit)
{
    float epsilon = 0.0001f;
    
    // ray direction = - view direction
    float3 reflectDirection = -2 * dot(ray.direction, hit.normal) * hit.normal + ray.direction;
    ray.direction = reflectDirection;
    // avoid acne: 
    ray.origin = hit.position + epsilon * hit.normal;
    ray.energy = ray.energy * hit.material.specular;
    return;
        
}

// Refracts the given ray from the given hit point
void refractRay(inout Ray ray, RayHit hit)
{
    // Your implementation
}

// Samples the _SkyboxTexture at a given direction vector
float3 sampleSkybox(float3 direction)
{
    float theta = acos(direction.y) / -PI;
    float phi = atan2(direction.x, -direction.z) / -PI * 0.5f;
    return _SkyboxTexture.SampleLevel(sampler_SkyboxTexture, float2(phi, theta), 0).xyz;
}
