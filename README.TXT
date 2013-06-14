The code contained in the Github repository SlowThought/SlowFlight (the current main application is zfoil.rkt; future references to zfoil refer to all software within the repository) is distributed under the GNU Lesser General Public License
(LGPL). This means that you can link zfoil into proprietary applications, provided you follow the rules stated in the LGPL. You can
also modify zfoil; if you distribute a modified version, you must distribute it under the terms of the LGPL, which in particular means
that you must publically release the source code for the modified software.

"zfoil.rkt" is an airfoil analysis tool. It predicts performance, given a geometry (as opposed to airfoil "design" tools, which propose
a geometry given a desired performance (coefficient, velocity distribution, whatever).

The non-linear nature of the equations governing fluid flow mean that zfoil, as any other computer code, can only estimate the true performance
of an airfoil. The user of zfoil assumes all risk. zfoil is provided without compensating consideration, "as is". The author hopes you don't die.
