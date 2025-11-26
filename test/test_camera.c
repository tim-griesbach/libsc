/*
  This file is part of the SC Library.
  The SC Library provides support for parallel scientific applications.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors

  The SC Library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  The SC Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with the SC Library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
*/

#include <sc_camera.h>

#define SC_CAMERA_TEST_EPS 1e-6

static void
check_difference(sc_camera_coords_t *expected,
                 sc_camera_coords_t *actual,
                 size_t length,
                 const char *msg)
{
  size_t i;
  for (i = 0; i < length; ++i)
  {
    SC_CHECK_ABORT(fabs(expected[i] - actual[i]) < SC_CAMERA_TEST_EPS,
                   msg);
  }
}

static void
check_difference_homogeneous( sc_camera_vec4_t expected, 
                              sc_camera_vec4_t actual, const char *msg)
{
  sc_camera_vec3_t scaled_expected, scaled_actual;

  scaled_expected[0] = expected[0] * actual[3];
  scaled_expected[1] = expected[1] * actual[3];
  scaled_expected[2] = expected[2] * actual[3];
  
  scaled_actual[0] = actual[0] * expected[3];
  scaled_actual[1] = actual[1] * expected[3];
  scaled_actual[2] = actual[2] * expected[3];

  check_difference(scaled_expected, scaled_actual, 3, msg);
}


static void 
quat_conjugate_transform(sc_camera_vec4_t out,
                                        const sc_camera_vec4_t q,
                                        const sc_camera_vec4_t p)
{
    // Unpack q
    const sc_camera_coords_t qx = q[0];
    const sc_camera_coords_t qy = q[1];
    const sc_camera_coords_t qz = q[2];
    const sc_camera_coords_t qw = q[3];

    // --- Step 1: temp = q * p ---
    const sc_camera_coords_t t_x = qw * p[0] + qx * p[3] + qy * p[2] - qz * p[1];
    const sc_camera_coords_t t_y = qw * p[1] - qx * p[2] + qy * p[3] + qz * p[0];
    const sc_camera_coords_t t_z = qw * p[2] + qx * p[1] - qy * p[0] + qz * p[3];
    const sc_camera_coords_t t_w = qw * p[3] - qx * p[0] - qy * p[1] - qz * p[2];

    // --- Step 2: q_inv = conjugate(q) = (-x, -y, -z, w) ---
    // out = temp * q_inv
    out[0] = t_w * (-qx) + t_x * qw + t_y * (-qz) - t_z * (-qy);
    out[1] = t_w * (-qy) - t_x * (-qz) + t_y * qw + t_z * (-qx);
    out[2] = t_w * (-qz) + t_x * (-qy) - t_y * (-qx) + t_z * qw;
    out[3] = t_w * qw - t_x * (-qx) - t_y * (-qy) - t_z * (-qz);
}

static void 
test_yaw_pitch_roll(sc_camera_t *camera)
{ 
  /* As in the sc_camera documentation defined the rotation quaternion (q*...*q^-1)
    is the rotation needed to move world points such that form the camera 
    perspective the camera rotates the world  */

  sc_camera_vec4_t p = {1., 1., 1., 0.};
  sc_camera_vec3_t expected;

  /* sets rotation to (0,0,0,1) (identity rotation) */
  sc_camera_init(camera);
  
  sc_camera_yaw(camera, M_PI/2.0);
  quat_conjugate_transform(p, camera->rotation, p);

  expected[0] = -1.; expected[1] = 1.; expected[2] = 1.;
  check_difference(expected, p, 3, "Yaw test failed.");

  p[0] = 1.; p[1] = 1.; p[2] = 1.; p[3] = 0.;

  sc_camera_init(camera);
  sc_camera_pitch(camera, M_PI/2.0);
  quat_conjugate_transform(p, camera->rotation, p);

  expected[0] = 1.; expected[1] = 1.; expected[2] = -1.;
  check_difference(expected, p, 3, "Pitch test failed.");

  p[0] = 1.; p[1] = 1.; p[2] = 1.; p[3] = 0.;

  sc_camera_init(camera);
  sc_camera_roll(camera, M_PI/2.0);
  quat_conjugate_transform(p, camera->rotation, p);

  expected[0] = 1.; expected[1] = -1.; expected[2] = 1.;
  check_difference(expected, p, 3, "Roll test failed.");
}

static void
test_look_at(sc_camera_t *camera)
{
  sc_camera_vec3_t eye = {1., 2., 3.};
  sc_camera_vec3_t center = {1., 0., -1.};
  sc_camera_vec3_t up = {0., 1., 0.};
  sc_camera_vec4_t p, expected;

  sc_camera_look_at(camera, eye, center, up);

  check_difference((sc_camera_vec3_t){1., 2., 3.}, camera->position, 3,
                   "Look at position test failed.");

  /* p = center - eye */
  p[0] = 0.; p[1] = -2.; p[2] = -4.; p[3] = 0.;
  quat_conjugate_transform(p, camera->rotation, p);

  expected[0] = 0.; expected[1] = 0.; expected[2] = -sqrt(20.); expected[3] = 0.;
  check_difference(expected, p, 4,
                   "Look at test failed.");

  /* p = up */
  p[0] = 0.; p[1] = 1.; p[2] = 0.; p[3] = 0.;
  quat_conjugate_transform(p, camera->rotation, p);

  expected[0] = 0.;
  check_difference(&p[0], expected, 1,
                   "Look at up x test failed.");
}

void
test_view_transform_trivial(sc_camera_t *camera)
{
  sc_camera_vec3_t in0 = {0., 0., 1.};
  sc_camera_vec3_t in1 = {2., 3., 5.};
  sc_camera_vec3_t expected;
  sc_array_t *points_in, *points_out;

  points_in  = sc_array_new_count(sizeof(sc_camera_vec3_t), 3);
  points_out = sc_array_new(sizeof(sc_camera_vec3_t));

  memcpy(sc_array_index(points_in, 0), &in0, sizeof(sc_camera_vec3_t));
  memcpy(sc_array_index(points_in, 1), &in1, sizeof(sc_camera_vec3_t));

  sc_camera_init(camera);

  sc_camera_view_transform(camera, points_in, points_out);

  expected[0] = 0.; expected[1] = 0.; expected[2] = 0.;
  check_difference (expected, sc_array_index(points_out, 0), 3, "View transform test failed.");

  expected[0] = 2.; expected[1] = 3.; expected[2] = 4.;
  check_difference (expected, sc_array_index(points_out, 1), 3, "View transform test failed.");

  sc_array_destroy(points_in);
  sc_array_destroy(points_out);
}

void 
test_view_transform(sc_camera_t *camera)
{
  sc_array_t *points_in, *points_out;
  sc_camera_vec3_t eye = {0.0, 0.0, 5.0};
  sc_camera_vec3_t center = {3.0,4.0,5.0};
  sc_camera_vec3_t up = {3.0, 4.0, 7.0};
  sc_camera_vec3_t expected;

  points_in  = sc_array_new_count(sizeof(sc_camera_vec3_t), 3);
  points_out = sc_array_new(sizeof(sc_camera_vec3_t));

  memcpy(sc_array_index(points_in, 0), &eye,    sizeof(sc_camera_vec3_t));
  memcpy(sc_array_index(points_in, 1), &center, sizeof(sc_camera_vec3_t));
  memcpy(sc_array_index(points_in, 2), &up,     sizeof(sc_camera_vec3_t));

  sc_camera_look_at(camera, eye, center, up);

  sc_camera_view_transform(camera, points_in, points_out);

  expected[0] = 0.; expected[1] = 0.; expected[2] = 0.; 
  check_difference(expected, sc_array_index(points_out, 0),
     3, "View transform test failed.");

  expected[0] = 0.; expected[1] = 0.; expected[2] = -5.; 
  check_difference(expected, sc_array_index(points_out, 1),
     3, "View transform test failed.");

  expected[0] = 0.; expected[1] = 2.; expected[2] = -5.; 
  check_difference(expected, sc_array_index(points_out, 2),
     3, "View transform test failed.");

  sc_array_destroy(points_in);
  sc_array_destroy(points_out);
}

void
test_get_view_mat(sc_camera_t *camera)
{
  sc_camera_mat4x4_t view_matrix;
  sc_camera_vec3_t eye = {-1., 4., 3.};
  sc_camera_vec3_t center = {2., 0., 1.};
  sc_camera_vec3_t up = {1., 1., 0.};

  /* The expected view matrix was calculated with GLM (OpenGl mathematics) library. */
  /* see : glm::lookAt<glm::f64> (The GLM library uses the same conventions for the view transform.) */
  sc_camera_mat4x4_t expected = {
    0.264906471413, 0.787070348709, -0.557086014531, 0.000000000000,
    -0.264906471413, 0.614898709929, 0.742781352708, 0.000000000000,
    0.927172649946, -0.049191896794, 0.371390676354, 0.000000000000,
    -1.456985592772, -1.524948800624, -4.642383454426, 1.000000000000
  };

  sc_camera_look_at(camera, eye, center, up);
  sc_camera_get_view_mat(camera, view_matrix);

  check_difference(expected,
                   view_matrix,
                   16,
                   "Get view matrix test failed.");
}

void 
test_projection_transform(sc_camera_t *camera)
{
  sc_array_t *points_in, *points_out;
  sc_camera_vec3_t in0 = {-0.01, 0.01, -0.01};
  sc_camera_vec3_t in1 = {100., -100., -100.};
  sc_camera_vec3_t in2 = {0., 0., -1.};
  sc_camera_vec4_t expected;
  sc_camera_coords_t *p;

  points_in  = sc_array_new_count(sizeof(sc_camera_vec3_t), 3);
  points_out = sc_array_new(sizeof(sc_camera_vec4_t));

  memcpy(sc_array_index(points_in, 0), &in0, sizeof(sc_camera_vec3_t));
  memcpy(sc_array_index(points_in, 1), &in1, sizeof(sc_camera_vec3_t));
  memcpy(sc_array_index(points_in, 2), &in2, sizeof(sc_camera_vec3_t));

  sc_camera_init(camera);
  sc_camera_projection_transform(camera, points_in, points_out);

  expected[0] = -1.; expected[1] = 1.; expected[2] = -1.; expected[3] = 1.;
  check_difference_homogeneous(expected, sc_array_index(points_out, 0), 
    "Projection transform test failed.");

  expected[0] = 1.; expected[1] = -1.; expected[2] = 1.; expected[3] = 1;
  check_difference_homogeneous(expected, sc_array_index(points_out, 1), 
    "Projection transform test failed.");

  p = sc_array_index(points_out, 2);
  SC_CHECK_ABORT(p[2] > -1. && p[2] < 1., "Projection transform test failed.");
  expected[0] = 0.; expected[1] = 0.; expected[2] = p[2]; expected[3] = p[3];
  check_difference_homogeneous(expected, p, "Projection transformation test failed.");

  sc_array_destroy(points_in);
  sc_array_destroy(points_out);
}

void
test_get_projection_mat(sc_camera_t *camera)
{
  sc_camera_mat4x4_t proj_matrix;

  /* The expected view matrix was calculated with GLM (OpenGl mathematics) library. */
  /* see : glm::perspective<glm::f64> (The GLM library uses slightly different conventions (FOV in y-direction).) */
  sc_camera_mat4x4_t expected = {
    1.732050807569, 0.000000000000, 0.000000000000, 0.000000000000,
    0.000000000000, 2.309401076759, 0.000000000000, 0.000000000000,
    0.000000000000, 0.000000000000, -1.002002002002, -1.000000000000,
    0.000000000000, 0.000000000000, -0.200200200200, 0.000000000000
  };

  sc_camera_aspect_ratio(camera, 8, 6);
  sc_camera_clipping_dist(camera, 0.1, 100.0);
  sc_camera_fov(camera, M_PI / 3.0);
  sc_camera_get_projection_mat(camera, proj_matrix);

  check_difference(expected,
                   proj_matrix,
                   16,
                   "Get projection matrix test failed.");
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_camera_t      *camera;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 0, 1, NULL, SC_LP_DEFAULT);

  camera = sc_camera_new();

  test_get_view_mat(camera);

  test_view_transform_trivial(camera);

  test_view_transform(camera);

  test_get_projection_mat(camera);

  test_projection_transform(camera);

  test_yaw_pitch_roll(camera);

  test_look_at(camera);

  sc_camera_destroy(camera);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}