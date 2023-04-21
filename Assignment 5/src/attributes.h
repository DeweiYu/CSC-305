#pragma once
///kkkkkk
#include <Eigen/Core>

class VertexAttributes
{
public:
    VertexAttributes(double x = 0, double y = 0, double z = 0, double w = 1)
    {
        position << x, y, z, w;
        normal << 0,0,0; //  initial value
    }

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(
        const VertexAttributes &a,
        const VertexAttributes &b,
        const VertexAttributes &c,
        const double alpha,
        const double beta,
        const double gamma)
    {
        VertexAttributes r;
        r.position = alpha * a.position + beta * b.position + gamma * c.position;
        r.normal = (a.normal + b.normal + c.normal).normalized();
        r.color = alpha * a.color + beta * b.color + gamma * c.color;
        return r;
    }

    Eigen::Vector3d normal;
    Eigen::Vector3d color;
    Eigen::Vector4d position;
   
};

class FragmentAttributes
{
public:
    FragmentAttributes(double r = 0, double g = 0, double b = 0, double a = 1)
    {
        color << r, g, b, a;
    }
    double depth = -2;
    Eigen::Vector4d color;
    Eigen::Vector4d position;
};

class FrameBufferAttributes
{
public:
    FrameBufferAttributes(uint8_t r = 0, uint8_t g = 0, uint8_t b = 0, uint8_t a = 255)
    {
        color << r, g, b, a;
    }
    double depth = -2;
    Eigen::Matrix<uint8_t, 4, 1> color;
};

class UniformAttributes
{
public:
    Eigen::Vector4d color;
    Eigen::Matrix4d view_trafos;
    bool pv = false;
};