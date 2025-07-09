#Building scifor
#Errors
set -e
set -v

cd scifor
pwd

mkdir -p build
cd build
pwd

echo "cmake .."
cmake ..

echo "make -j"
make -j

echo "make install"
make install

CONF_PATH=$(grep CONFIG_FILE_PATH CMakeCache.txt |cut -d "=" -f2)
echo "source $CONF_PATH/scifor_config_user.sh" >> ~/.scifor_config_user
echo -e "\e[32m scifor compiled and sourced \e[0m"
