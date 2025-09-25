# fluor-fos
Open-source software (Fluor-FOS) for fluorescent nanocomposite multi-layer optical modeling. The program can handel infinte number of layers, scattering particlles, and fluorescne tparticles, all with distnct material and optical properteis. This progran is equpped with user-freindly interface that allows the user to insert the particles and mediums and then place then in multi-layer structure with adjustable struicutral parameters. The input file (input.txt) can be inserted once the main progrma file is excuted. The program is in python-language and can be downlaoded to be run via any operating system, along with excutable file. 
You can refer to our article for more details (). 

The excutable file can be found in the resciporty .

The main functionalities our program offer can be summarized as follow:
  1) Mie theory calcualtion for spherical aprticle inclusions. For non-spherical inclusion, the user can input already calcualted spectral attenutaition properties (lambda, mu_s, mu_a, and g) by other numerical method such as COMSOL. Other cases like dye inclusions has been discussed in the article.
  2) The user can insert multiple fluorescent particle inclusions in multi-layer with various light incednet sources for creating new colors or spectral shifting purposes for 1D systems, such strategy can be found in radiaitve cooling and energy harvesting applications. The software has the potential to be expanded to 2D and 3D systems as future develpments 
  3) The user can insert uncertain optical properties by providing the standard deviation along with optical porpertei and how many runs desired for smooth intervals. The uncertaintly analysis is based on a derived equation of normal distrbution coupled with random number via the fundamnetal law of Monte Carlo ().
  4) The software has been biult utilizing parallel processing for photon launching using Numba pacakge, despite the difficulties with 2D tensor reflectance, transmittance, and absorptance.
  5) The metalic or polymeric substrates can be included.



