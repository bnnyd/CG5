342476157 bennydv
204880421 ellavidan

Q4.8:
Given the plane normal, check plane's orientation (XY, YZ, ZX). Define the two coordinates as u and v. 
Check if each coordinate (u and v) is closer to its floor() or to its floor() + 1. 
If u and v agree on the side, define the material m1, else m2. The agreement is deduced by using the  xor logical operation (it's True when the coordinates disagree)

Q.5:
In the infinite cylinder formula, replace x and z with a point on the current ray's x and z, i.e. (orig + dir * t).x or .z
If there is no valid intersection between the ray and the infinite cylinder there will be no intersection with the finite cylinder.
If there is an intersection, check if intersection y coordinate is in a valid the height. If yes, return the current ray hit. Else, check if the ray intersects the caps.