// This file was generated by ODB, object-relational mapping (ORM)
// compiler for C++.
//

#ifndef IMAGE_ENTITY_ODB_HXX
#define IMAGE_ENTITY_ODB_HXX

#include <odb/version.hxx>

#if (ODB_VERSION != 10100UL)
#error ODB runtime version mismatch
#endif

#include <odb/pre.hxx>

// Begin prologue.
//
//
// End prologue.

#include "ImageEntity.hxx"

#include <memory>
#include <cstddef>

#include <odb/core.hxx>
#include <odb/traits.hxx>
#include <odb/pointer-traits.hxx>
#include <odb/container-traits.hxx>
#include <odb/result.hxx>

#include <odb/mysql/version.hxx>
#include <odb/mysql/forward.hxx>
#include <odb/mysql/mysql-types.hxx>
#include <odb/mysql/query.hxx>

#include <odb/details/buffer.hxx>

namespace odb
{
  // ImageEntity
  //
  template <>
  class access::object_traits< ::ImageEntity >
  {
    public:
    typedef ::ImageEntity object_type;
    typedef ::ImageEntity* pointer_type;
    typedef long unsigned int id_type;

    struct image_type
    {
      // id_
      //
      unsigned long long id_value;
      my_bool id_null;

      // name_
      //
      details::buffer name_value;
      unsigned long name_size;
      my_bool name_null;

      // study_
      //
      details::buffer study_value;
      unsigned long study_size;
      my_bool study_null;

      // measureid_
      //
      int measureid_value;
      my_bool measureid_null;

      // localid_
      //
      int localid_value;
      my_bool localid_null;

      // value_
      //
      double value_value;
      my_bool value_null;

      // lambda_
      //
      double lambda_value;
      my_bool lambda_null;

      // expid_
      //
      int expid_value;
      my_bool expid_null;

      // generation_
      //
      int generation_value;
      my_bool generation_null;

      // slice_
      //
      int slice_value;
      my_bool slice_null;

      std::size_t version;
    };

    struct id_image_type
    {
      unsigned long long id_value;
      my_bool id_null;

      std::size_t version;
    };

    typedef mysql::query query_base_type;

    struct query_type: query_base_type
    {
      // id
      //
      static const mysql::query_column<
        mysql::value_traits< long unsigned int, unsigned long long, mysql::id_ulonglong >::query_type,
        mysql::id_ulonglong>
      id;

      // name
      //
      static const mysql::query_column<
        mysql::value_traits< ::std::string, details::buffer, mysql::id_string >::query_type,
        mysql::id_string>
      name;

      // study
      //
      static const mysql::query_column<
        mysql::value_traits< ::std::string, details::buffer, mysql::id_string >::query_type,
        mysql::id_string>
      study;

      // measureid
      //
      static const mysql::query_column<
        mysql::value_traits< int, int, mysql::id_long >::query_type,
        mysql::id_long>
      measureid;

      // localid
      //
      static const mysql::query_column<
        mysql::value_traits< int, int, mysql::id_long >::query_type,
        mysql::id_long>
      localid;

      // value
      //
      static const mysql::query_column<
        mysql::value_traits< double, double, mysql::id_double >::query_type,
        mysql::id_double>
      value;

      // lambda
      //
      static const mysql::query_column<
        mysql::value_traits< double, double, mysql::id_double >::query_type,
        mysql::id_double>
      lambda;

      // expid
      //
      static const mysql::query_column<
        mysql::value_traits< int, int, mysql::id_long >::query_type,
        mysql::id_long>
      expid;

      // generation
      //
      static const mysql::query_column<
        mysql::value_traits< int, int, mysql::id_long >::query_type,
        mysql::id_long>
      generation;

      // slice
      //
      static const mysql::query_column<
        mysql::value_traits< int, int, mysql::id_long >::query_type,
        mysql::id_long>
      slice;

      query_type ();
      query_type (const std::string&);
      query_type (const query_base_type&);
    };

    static const std::size_t in_column_count = 10UL;
    static const std::size_t out_column_count = 10UL;

    static const char* const persist_statement;
    static const char* const find_statement;
    static const char* const update_statement;
    static const char* const erase_statement;
    static const char* const query_clause;

    struct container_statement_cache_type;

    static id_type
    id (const object_type&);

    static id_type
    id (const image_type&);

    static void
    grow (image_type&, my_bool*);

    static void
    bind (MYSQL_BIND*, image_type&, bool);

    static void
    bind (MYSQL_BIND*, id_image_type&);

    static void
    init (image_type&, const object_type&);

    static void
    init (object_type&, const image_type&, database&);

    static void
    init (id_image_type&, const id_type&);

    static void
    persist (database&, object_type&);

    static void
    update (database&, const object_type&);

    static void
    erase (database&, const id_type&);

    static pointer_type
    find (database&, const id_type&);

    static bool
    find (database&, const id_type&, object_type&);

    template<typename T>
    static result<T>
    query (database&, const query_type&);

    public:
    static bool
    find_ (mysql::object_statements< object_type >&, const id_type&);

    static void
    load_ (mysql::object_statements< object_type >&, object_type&);

    static void
    query_ (database&,
            const query_type&,
            mysql::object_statements< object_type >&,
            details::shared_ptr< mysql::select_statement >&);
  };
}

#include "ImageEntity-odb.ixx"

// Begin epilogue.
//
//
// End epilogue.

#include <odb/post.hxx>

#endif // IMAGE_ENTITY_ODB_HXX