
#include <car.h>

const char         *car_type = "car";

static const char  *
get_type_fn (sc_object_t * o)
{
  SC_LDEBUG ("car get_type\n");

  return car_type;
}

static void
initialize_fn (sc_object_t * o)
{
  Car                *car = car_get_data (o);

  SC_LDEBUG ("car initialize\n");

  car->speed = 0;
  car->wheelsize = 17;
}

static void
write_fn (sc_object_t * o, FILE * out)
{
  Car                *car = car_get_data (o);

  fprintf (out, "Car (wheel size %f) speeds at %f km/h\n",
           car->wheelsize, car->speed);
}

static float
wheelsize_fn (sc_object_t * o)
{
  Car                *car = car_get_data (o);

  SC_LDEBUG ("car wheelsize\n");

  return car->wheelsize;
}

sc_object_t        *
car_klass_new (sc_object_t * d)
{
  bool                a1, a2, a3, a4;
  sc_object_t        *o;

  SC_ASSERT (d != NULL);

  o = sc_object_alloc ();
  sc_object_delegate_push (o, d);

  a1 = sc_object_method_register (o, (sc_object_method_t) sc_object_get_type,
                                  (sc_object_method_t) get_type_fn);
  a2 =
    sc_object_method_register (o, (sc_object_method_t) sc_object_initialize,
                               (sc_object_method_t) initialize_fn);
  a3 =
    sc_object_method_register (o, (sc_object_method_t) sc_object_write,
                               (sc_object_method_t) write_fn);
  a4 =
    sc_object_method_register (o, (sc_object_method_t) car_wheelsize,
                               (sc_object_method_t) wheelsize_fn);
  SC_ASSERT (a1 && a2 && a3 && a4);

  sc_object_initialize (o);

  return o;
}

Car                *
car_get_data (sc_object_t * o)
{
  SC_ASSERT (sc_object_is_type (o, car_type));

  return (Car *) sc_object_get_data (o, (sc_object_method_t) car_get_data,
                                     sizeof (Car));
}

float
car_wheelsize (sc_object_t * o)
{
  sc_object_method_t  oinmi;

  oinmi = sc_object_delegate_lookup (o, (sc_object_method_t) car_wheelsize);
  SC_ASSERT (oinmi != NULL);

  return ((float (*)(sc_object_t *)) oinmi) (o);
}