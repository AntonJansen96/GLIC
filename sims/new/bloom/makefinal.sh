#!/bin/bash

for i in 1 2 3 4 5 6; do convert fig${i}_left.png fig${i}_right.png +append fig${i}_comb.png; done
convert fig[123]_comb.png -append A.png
convert fig[456]_comb.png -append B.png
convert A.png B.png +append final.png
rm A.png B.png
