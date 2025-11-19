#
# Update module list page
# 

module load Anaconda3/2023.09-0

echo "Building new module list..."

python modpage.py

echo "modules.rst file built - check for errors then run push.sh"


