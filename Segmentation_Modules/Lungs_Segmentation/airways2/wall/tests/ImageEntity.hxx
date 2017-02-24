#ifndef IMAGEENTITY_HXX
#define IMAGEENTITY_HXX
#include <string>

#include <odb/core.hxx>     // (1)

#pragma db object           // (2)
class ImageEntity
{
public:

	ImageEntity (const std::string& name, const std::string& study, int mid, int lid, double value, double lambda, int expID, int generation)
	: name_ (name), study_ (study), measureid_ (mid), localid_ (lid), value_ (value), lambda_ (lambda), expid_ (expID), generation_(generation)
	{
	}
	ImageEntity (const std::string& name, const std::string& study, int mid, int lid, double value, double lambda, int expID, int generation, int slice)
	: name_ (name), study_ (study), measureid_ (mid), localid_ (lid), value_ (value), lambda_ (lambda), expid_ (expID), generation_(generation), slice_ (slice)
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
  ImageEntity () {}              // (3)

  friend class odb::access; // (4)

  #pragma db id auto        // (5)
  unsigned long id_;        // (5)

  #pragma db key_type("VARCHAR(255) NOT NULL")
  std::string name_;
  #pragma db key_type("VARCHAR(255) NOT NULL")
  std::string study_;
  #pragma db key_type("INT UNSIGNED NOT NULL")
  int measureid_;
  #pragma db key_type("INT UNSIGNED NOT NULL")
  int localid_;

  double value_;
  #pragma db key_type("DOUBLE NOT NULL")
  double lambda_;
  #pragma db key_type("INT UNSIGNED NOT NULL")
  int expid_;
  #pragma db key_type("INT UNSIGNED NOT NULL")
  int generation_;
  #pragma db key_type("INT UNSIGNED")
  int slice_;

};
#endif
