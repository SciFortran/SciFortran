#!/bin/sh
# Install the artefacts of a SciFortran build into a Debian staging tree.
#
# Upstream's "make install" is not usable for packaging (it writes into $HOME,
# uses a non-FHS layout and ships environment-modules files), so the three
# things that matter are installed by hand:
#
#   libscifor.a       -> usr/lib/<multiarch>/
#   *.mod             -> usr/lib/<multiarch>/fortran/gfortran-mod-<N>/scifor/
#   scifor.pc         -> usr/lib/<multiarch>/pkgconfig/   (written here)
#
# <N> is the gfortran module-format version read from the .mod files; the same
# value is published as the ${fortran:Depends} substvar (gfortran-mod-<N>), so
# the package cannot be installed next to a compiler that cannot read its
# modules.
#
# Usage: scifor-install.sh BUILDDIR STAGEDIR MULTIARCH UPSTREAM_VERSION SUBSTVARS

set -eu

if [ "$#" -ne 5 ]; then
    echo "usage: $0 BUILDDIR STAGEDIR MULTIARCH UPSTREAM_VERSION SUBSTVARS" >&2
    exit 2
fi
build=$1
stage=$2
multiarch=$3
version=$4
substvars=$5

die() { echo "scifor-install: error: $*" >&2; exit 1; }

lib=$build/libscifor.a
incdir=$build/include
genpc=$build/etc/scifor.pc

[ -f "$lib" ]    || die "$lib not found (did the build succeed?)"
[ -d "$incdir" ] || die "$incdir not found"
[ -f "$genpc" ]  || die "$genpc not found"
ls "$incdir"/*.mod >/dev/null 2>&1 || die "no .mod files in $incdir"

# gfortran module-format version, e.g. "GFORTRAN module version '15' created ..."
first_mod=$(ls "$incdir"/*.mod | head -n 1)
modver=$(zcat -f "$first_mod" | head -n 1 \
         | sed -n "s/^GFORTRAN module version '\([0-9][0-9]*\)'.*/\1/p")
[ -n "$modver" ] || die "cannot determine the gfortran module version from $first_mod"

libdir=usr/lib/$multiarch
fmoddir=$libdir/fortran/gfortran-mod-$modver/scifor
pcdir=$libdir/pkgconfig

install -d "$stage/$libdir" "$stage/$fmoddir" "$stage/$pcdir"
install -m 0644 "$lib" "$stage/$libdir/libscifor.a"
install -m 0644 "$incdir"/*.mod "$stage/$fmoddir/"

# Link line: take what upstream's CMake worked out for the dependencies
# (ScaLAPACK, LAPACK, BLAS) but turn absolute paths into plain -l flags and
# drop the system library directories, so that the .pc file does not depend on
# the build machine's layout and honours the BLAS/LAPACK alternatives.
extra=
for tok in $(sed -n 's/^Libs:[[:space:]]*//p' "$genpc" | tr ';' ' '); do
    case $tok in
        -lscifor|-L'${libdir}'|-L/usr/lib|-L/usr/lib/"$multiarch") ;;
        /*/lib*.so|/*/lib*.a)
            name=${tok##*/}; name=${name#lib}; name=${name%.*}
            extra="$extra -l$name" ;;
        *) extra="$extra $tok" ;;
    esac
done

# Fail loudly if a dependency silently dropped out of the build.
for dep in scalapack lapack blas; do
    case $extra in
        *-l$dep*) ;;
        *) die "-l$dep is missing from the generated link line:$extra" ;;
    esac
done

cat > "$stage/$pcdir/scifor.pc" <<PC
prefix=/usr
libdir=\${prefix}/lib/$multiarch
fmoddir=\${libdir}/fortran/gfortran-mod-$modver/scifor

Name: scifor
Description: SciFortran, a Fortran library for mathematics, science and engineering
URL: https://github.com/SciFortran/SciFortran
Version: $version
Cflags: -I\${fmoddir}
Libs: -L\${libdir} -lscifor$extra
PC
chmod 0644 "$stage/$pcdir/scifor.pc"

# ${fortran:Depends}, read by dh_gencontrol for libscifor-dev.
printf 'fortran:Depends=gfortran-mod-%s\n' "$modver" > "$substvars"
