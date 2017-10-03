#include<libgeodecomp/config.h>
#ifdef LIBGEODECOMP_WITH_BOOST_SERIALIZATION
#ifndef LIBGEODECOMP_BOOSTSERIALIZATION_H
#define LIBGEODECOMP_BOOSTSERIALIZATION_H



namespace LibGeoDecomp {
class BoostSerialization
{
public:
<<<<<<< HEAD
=======
    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::Chronometer& object, const unsigned /*version*/)
    {
        archive & object.totalTimes;
    }

    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::Color& object, const unsigned /*version*/)
    {
        archive & object.rgb;
    }

    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::Coord<1 >& object, const unsigned /*version*/)
    {
        archive & object.c;
    }

    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::Coord<2 >& object, const unsigned /*version*/)
    {
        archive & object.c;
    }

    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::Coord<3 >& object, const unsigned /*version*/)
    {
        archive & object.c;
    }

    template<typename ARCHIVE, int DIM>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::CoordBox<DIM>& object, const unsigned /*version*/)
    {
        archive & object.dimensions;
        archive & object.origin;
    }

    template<typename ARCHIVE, typename T, int SIZE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::FixedArray<T, SIZE>& object, const unsigned /*version*/)
    {
        archive & object.elements;
        archive & object.store;
    }

    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::FloatCoord<1 >& object, const unsigned /*version*/)
    {
        archive & object.c;
    }

    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::FloatCoord<2 >& object, const unsigned /*version*/)
    {
        archive & object.c;
    }

    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::FloatCoord<3 >& object, const unsigned /*version*/)
    {
        archive & object.c;
    }

    template<typename ARCHIVE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::NonPoDTestCell& object, const unsigned /*version*/)
    {
        archive & object.coord;
        archive & object.cycleCounter;
        archive & object.missingNeighbors;
        archive & object.seenNeighbors;
        archive & object.simSpace;
    }

    template<typename ARCHIVE, typename VALUE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::Palette<VALUE>& object, const unsigned /*version*/)
    {
        archive & object.colors;
    }

    template<typename ARCHIVE, typename VALUE>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::QuickPalette<VALUE>& object, const unsigned /*version*/)
    {
        archive & object.mark0;
        archive & object.mark1;
        archive & object.mark2;
        archive & object.mark3;
        archive & object.mark4;
        archive & object.mult;
    }

    template<typename ARCHIVE, int DIMENSIONS>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::Region<DIMENSIONS>& object, const unsigned /*version*/)
    {
        archive & object.geometryCacheTainted;
        archive & object.indices;
        archive & object.myBoundingBox;
        archive & object.mySize;
    }

    template<typename ARCHIVE, int DIM>
    inline
    static void serialize(ARCHIVE& archive, LibGeoDecomp::Streak<DIM>& object, const unsigned /*version*/)
    {
        archive & object.endX;
        archive & object.origin;
    }
>>>>>>> ba45a35aa1efd409063055b54a253f5a65944985

};
}

namespace boost {
namespace serialization {

using namespace LibGeoDecomp;

<<<<<<< HEAD
=======
template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::Chronometer& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::Color& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::Coord<1 >& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::Coord<2 >& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::Coord<3 >& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE, int DIM>
void serialize(ARCHIVE& archive, LibGeoDecomp::CoordBox<DIM>& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE, typename T, int SIZE>
void serialize(ARCHIVE& archive, LibGeoDecomp::FixedArray<T, SIZE>& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::FloatCoord<1 >& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::FloatCoord<2 >& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::FloatCoord<3 >& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE>
void serialize(ARCHIVE& archive, LibGeoDecomp::NonPoDTestCell& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE, typename VALUE>
void serialize(ARCHIVE& archive, LibGeoDecomp::Palette<VALUE>& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE, typename VALUE>
void serialize(ARCHIVE& archive, LibGeoDecomp::QuickPalette<VALUE>& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE, int DIMENSIONS>
void serialize(ARCHIVE& archive, LibGeoDecomp::Region<DIMENSIONS>& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

template<class ARCHIVE, int DIM>
void serialize(ARCHIVE& archive, LibGeoDecomp::Streak<DIM>& object, const unsigned version)
{
    BoostSerialization::serialize(archive, object, version);
}

>>>>>>> ba45a35aa1efd409063055b54a253f5a65944985

}
}

#endif

#endif
