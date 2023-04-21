// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"

#include <gif.h>
#include <fstream>

#include <Eigen/Geometry>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

using namespace std;
using namespace Eigen;

//Image height
const int H = 480;

//Camera settings
const double near_plane = 1.5; //AKA focal length
const double far_plane = near_plane * 100;
const double field_of_view = 0.7854; //45 degrees
const double aspect_ratio = 1.5;
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 3);
const Vector3d camera_gaze(0, 0, -1);
const Vector3d camera_top(0, 1, 0);

//Object
const std::string data_dir = DATA_DIR;
const std::string mesh_filename(data_dir + "bunny.off");
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)

//Material for the object
const Vector3d obj_diffuse_color(0.5, 0.5, 0.5);
const Vector3d obj_specular_color(0.2, 0.2, 0.2);
const double obj_specular_exponent = 256.0;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector3d> light_colors;
//Ambient light
const Vector3d ambient_light(0.3, 0.3, 0.3);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    if (!in.good())
    {
        std::cerr << "Invalid file " << mesh_filename << std::endl;
        exit(1);
    }
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16);
}

void build_uniform(UniformAttributes &uniform)
{
    //TODO: setup uniform

    //TODO: setup camera, compute w, u, v
    const Vector3d w = -camera_gaze.normalized();
    const Vector3d u = camera_top.cross(w).normalized();
    const Vector3d v = w.cross(u);
    

    //TODO: compute the camera transformation
    Matrix4d Matrix_cam;

    Matrix_cam <<
        u(0), v(0), w(0), camera_position(0),u(1), v(1), w(1), camera_position(1),u(2), v(2), w(2), camera_position(2),
        0, 0, 0, 1;
    Matrix_cam = Matrix_cam.inverse().eval();;
    
    //TODO: setup projection matrix
    
    Matrix4d Morth;
    
    const Vector3d close = camera_position + near_plane * camera_gaze;
    const Vector3d furtherest_point = camera_position + far_plane * camera_gaze;
    
    double y = near_plane * tan(field_of_view / 2);
    double x = near_plane * tan(field_of_view / 2)* aspect_ratio;
    
    Vector3d close_corner = close - y * v - x * u;
    Vector3d far_corner = furtherest_point + y * v + x * u;
    
   double l, b, n, r, t, f;
    l = close_corner(0);
    b = close_corner(1);
    n = close_corner(2);
    
    r = far_corner(0);
    t = far_corner(1);
    f = far_corner(2);

    Morth <<
        2 / (r - l), 0, 0, -(r + l) / (r - l),
        0, 2 / (t - b), 0, -(t + b) / (t - b),
        0, 0, 2 / (n - f), -(n + f) / (n - f),
        0, 0, 0, 1;
    
    Matrix4d P;
    
    P <<
        n, 0, 0, 0,
        0, n, 0, 0,
        0, 0, n + f, -f * n,
        0, 0, 1, 0;
    
    if (is_perspective)
    {
        //TODO setup prespective camera
        uniform.view_trafos = Morth * Matrix_cam; 
    }
    else
    {
        uniform.view_trafos = Morth * Matrix_cam;
    }
}

void simple_render(Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader

        VertexAttributes output;
        output.position = uniform.view_trafos * va.position;
        return output;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        //FragmentAttributes out(uniform.color(0), uniform.color(1), uniform.color(2));
        //out.depth = va.position(2);
        //return out;
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: fill the shader
        
 
        
        return FrameBufferAttributes(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255);
    };
    
    std::vector<VertexAttributes> vertex_attributes;
    //TODO: build the vertex attributes from vertices and facets

  
    for (int i = 0; i < facets.rows(); ++i) {
    // Get the vertex coordinates for the three vertices of the current facet and add them to vertex_attributes
        for (int j = 0; j < 3; ++j) {
        int vertex_index = facets(i, j);
        float x = vertices(vertex_index, 0);
        float y = vertices(vertex_index, 1);
        float z = vertices(vertex_index, 2);
        vertex_attributes.push_back(VertexAttributes(x, y, z));
    }
}

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
    
}
Matrix4d compute_rotation(const double alpha)
{
    //TODO: Compute the rotation matrix of angle alpha on the y axis around the object barycenter
// Create a rotation matrix around the y-axis with angle alpha
    Matrix4d rotation_matrix;
    double sin_alpha = sin(alpha);
    double cos_alpha = cos(alpha);
    
    rotation_matrix <<
     cos_alpha, 0, sin_alpha, 0,
    0, 1, 0, 0,
    -sin_alpha, 0, cos_alpha, 0,
    0, 0, 0, 1;


    return rotation_matrix;

}

void wireframe_render(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    Matrix4d trafo = compute_rotation(alpha);

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        VertexAttributes output;
        output.position = uniform.view_trafos * va.position;
        return output;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: fill the shader
        return FrameBufferAttributes(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255);
    };

    std::vector<VertexAttributes> vertex_attributes;

    //TODO: generate the vertex attributes for the edges and rasterize the lines
    //TODO: use the transformation matrix
    
    vertex_attributes.reserve(facets.rows() * 6); // pre-allocate memory
    for (int i = 0; i < facets.rows(); ++i) {
    const int v0 = facets(i, 0);
    const int v1 = facets(i, 1);
    const int v2 = facets(i, 2);
    vertex_attributes.emplace_back(vertices(v0, 0), vertices(v0, 1), vertices(v0, 2));
    vertex_attributes.emplace_back(vertices(v1, 0), vertices(v1, 1), vertices(v1, 2));
    vertex_attributes.emplace_back(vertices(v1, 0), vertices(v1, 1), vertices(v1, 2));
    vertex_attributes.emplace_back(vertices(v2, 0), vertices(v2, 1), vertices(v2, 2));
    vertex_attributes.emplace_back(vertices(v2, 0), vertices(v2, 1), vertices(v2, 2));
    vertex_attributes.emplace_back(vertices(v0, 0), vertices(v0, 1), vertices(v0, 2));
}

    rasterize_lines(program, uniform, vertex_attributes, 0.5, frameBuffer);

    std::vector<uint8_t> image(frameBuffer.size() * 3);
    const char *fileName = "wireflame.gif";
    int delay = 10;
    GifWriter g;
    GifBegin(&g, fileName, frameBuffer.rows(), frameBuffer.cols(), delay);

for (double t = 0; t < 15; t++) {
    frameBuffer.setConstant(FrameBufferAttributes());
    uniform.view_trafos *= trafo;
    rasterize_lines(program, uniform, vertex_attributes, 1, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
}

GifEnd(&g);
}

void get_shading_program(Program &program)
{
    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: transform the position and the normal
        //TODO: compute the correct lighting
        VertexAttributes output;
        output.position = uniform.view_trafos * va.position;
        output.normal = va.normal;
        output.color = va.color;

        return output;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        FragmentAttributes output;
        if (!uniform.pv) {
        Vector3d p(va.position(0), va.position(1), va.position(2));
        Vector3d N = va.normal;
        Vector3d ray_direction = p - camera_position;
        Vector3d lights_color;

    for (int i = 0; i < light_positions.size(); ++i) {
        const Vector3d& light_position = light_positions[i];
        const Vector3d& light_color = light_colors[i];

        const Vector3d Li = (light_position - p).normalized();

        Vector3d diff_color = obj_diffuse_color;

        const Vector3d diffuse = diff_color * std::max(Li.dot(N), 0.0);

        Vector3d v = -ray_direction.normalized();// notes from a3/a4
        Vector3d h = (Li + v).normalized();
        const Vector3d specular = obj_specular_color * pow(std::max(h.dot(N), 0.), obj_specular_exponent);

        const Vector3d D = light_position - p;
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }
    lights_color = lights_color + ambient_light;
    output.color << lights_color(0), lights_color(1), lights_color(2);
            }
        else {
        output.color << va.color(0), va.color(1), va.color(2);
        }

        output.position = va.position;
            return output;

    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: implement the depth check
        FrameBufferAttributes output;
        if (fa.position[2] > previous.depth) {
            output.color[0] = static_cast<uint8_t>(fa.color[0] * 255);
            output.color[1] = static_cast<uint8_t>(fa.color[1] * 255);
            output.color[2] = static_cast<uint8_t>(fa.color[2] * 255);
            output.color[3] = static_cast<uint8_t>(fa.color[3] * 255);
            output.depth = fa.position[2];
            return output;
        } else {
            return previous;
}

    };
}

void flat_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);
    Eigen::Matrix4d trafo = compute_rotation(alpha);

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: compute the normals
    //TODO: set material colors

    

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

void pv_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);

    Eigen::Matrix4d trafo = compute_rotation(alpha);

    //TODO: compute the vertex normals as vertex normal average

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: create vertex attributes
    //TODO: set material colors

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

int main(int argc, char *argv[])
{
    setup_scene();

    int W = H * aspect_ratio;
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer(W, H);
    vector<uint8_t> image;

    simple_render(frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("simple.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer_wire(W, H);
   
    wireframe_render(0, frameBuffer_wire);
    framebuffer_to_uint8(frameBuffer_wire, image);
    stbi_write_png("wireframe.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    flat_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("flat_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
   
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer_new(W, H);
    pv_shading(0, frameBuffer_new);
    framebuffer_to_uint8(frameBuffer_new, image);
    stbi_write_png("pv_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    //TODO: add the animation
     flat_shading(1, frameBuffer);
     pv_shading(1, frameBuffer);
     wireframe_render(1, frameBuffer);
    return 0;
}
