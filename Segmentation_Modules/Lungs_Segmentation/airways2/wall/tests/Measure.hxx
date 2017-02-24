#ifndef MEASURE_HXX
#define MEASURE_HXX
#include <string>

#include <odb/core.hxx>     // (1)

#pragma db object           // (2)
class Measure
{
public:
	  Measure (const std::string& name)
	      : name_ (name)
	  {
	  }

	  const std::string&
	  name () const
	  {
	    return name_;
	  }

	  unsigned long
	  id ()
	  {
	    return id_;
	  }

	  void
	  name (const std::string& name)
	  {
	    name_ = name;
	  }

private:
  Measure () {}              // (3)

  friend class odb::access; // (4)

  #pragma db id auto        // (5)
  unsigned long id_;        // (5)

  std::string name_;
};
#endif // MEASURE_HXX
