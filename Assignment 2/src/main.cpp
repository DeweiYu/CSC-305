// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h 引用library
using namespace Eigen;

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask 

    const Vector3d camera_origin(0, 0, 3);// original of the camera
    const Vector3d camera_view_direction(0, 0, -1);// the original of the direction of the camera

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1); //image position
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);// the displacemet of X
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);//the displacement for y

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    
    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement; // determine the pixel of the pixel center piece

            // Prepare the ray
            
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction =  camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;
            
             
            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                //
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png

    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    
    const Vector3d pgram_origin(-0.5, -0.5, 0);  // original position
    const Vector3d pgram_u(0, 0.7, -10);   // x- directionn
    const Vector3d pgram_v(1, 0.4, 0);    // y-direction

    
    
    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
           
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // TODO: Check if the ray intersects with the parallelogram
            Matrix3d scale;// formula for ray-tracing arallelogram to check if it intersect
            scale <<pgram_u[0],pgram_v[0],-ray_direction[0],pgram_u[1],pgram_v[1],-ray_direction[1],
            pgram_u[2],pgram_v[2],-ray_direction[2];
            
            Vector3d y_position = ray_origin - pgram_origin;
            Vector3d radius = scale.inverse()*y_position;
        

            if (0<=radius[0] && radius[0]<=1 && 0<=radius[1] && radius[1]<=1)
            {
                // TODO: The ray hit the parallelogram, compute the exact intersection
                
                Vector3d ray_intersection = ray_origin + ray_direction * radius[2];

                // TODO: Compute normal at the intersection point
                Vector3d normal_vector = (-1)*pgram_u.cross(pgram_v);  // reverse the direction 
                Vector3d ray_normal = normal_vector.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800);
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = ray_origin - camera_origin;
            
            Matrix3d scale;// copy from te prallagram function above to check
            scale<<pgram_u[0],pgram_v[0],-ray_direction[0],pgram_u[1],pgram_v[1],-ray_direction[1],
                    pgram_u[2],pgram_v[2],-ray_direction[2];
            Vector3d y_position = ray_origin - pgram_origin;
            
            Vector3d radius = scale.inverse()*y_position;
           
           
            // TODO: Check if the ray intersects with the parallelogram
            if (0<=radius[0] && radius[0]<=1 && 0<=radius[1] && radius[1]<=1)
            {
                // TODO: The ray hit the parallelogram, compute the exact intersection
                
                Vector3d ray_intersection = ray_origin + ray_direction * radius[2];

                // TODO: Compute normal at the intersection point
                //only one normal for the parallelogram
                Vector3d normal_vector = (-1)*pgram_u.cross(pgram_v);
                Vector3d ray_normal = normal_vector.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    // MatrixXd C= MatrixXd::Zero(800, 800); // Store the color
    MatrixXd R = MatrixXd::Zero(800, 800); 
    MatrixXd G = MatrixXd::Zero(800, 800);
    MatrixXd B = MatrixXd::Zero(800, 800);
    //MatrixXd C = MatrixXd::Zero(800, 800); 
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < R.cols(); ++i)
    {
        for (unsigned j = 0; j < R.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            
            const double sphere_radius = 0.9; 

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // TODO: Add shading parameter here
                
                // declaration of l v n, light direction,  view direction and the cross line made 90 degree
                Vector3d l = (light_position - ray_intersection).normalized();
                Vector3d v = -ray_direction.normalized();
                Vector3d n = (l+v).normalized();
                
               
                const double diffuse = std::max(l.dot(ray_normal),0.0);
                const double specular = pow(std::max(n.dot(ray_normal),0.0), specular_exponent);

                // Simple diffuse model
                R(i, j) = ambient + diffuse_color[0]*diffuse + specular_color[0]*specular;
                G(i, j) = ambient + diffuse_color[1]*diffuse + specular_color[1]*specular;
                B(i, j) = ambient + diffuse_color[2]*diffuse + specular_color[2]*specular;
                // implemented the red green blue model baes of the formular ambinet + diffuser + specular;

                // Clamp to zero 
                R(i, j) = std::max(R(i, j), 0.);
                G(i, j) = std::max(G(i, j), 0.);
                B(i, j) = std::max(B(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}
