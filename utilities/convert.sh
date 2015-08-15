#!/bin/bash

diff A.out A.out.reference
diff apply.out apply.out.reference
diff LU.out LU.out.reference


./removeParens.sh matrix.out
./deal2MTX LU.out LU.MTX
./deal2MTX matrix.out matrix.MTX
./deal2MTX A.out A.MTX
