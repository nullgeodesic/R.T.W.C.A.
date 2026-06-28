# About
RTWCA is a modular general relativistic ray tracer written in Julia for producing true color images in semi-arbitrary spacetimes. It runs on CPUs and GPUs (for now, only Nvidia's).

## Features:
* Integrators read from source files specifying everything specific about the scene to be rendered, then compile the code, so that all the scenes can be rendered with the same integration code, and new ones can be coded qwickly, with the others as a template.
* Uses adaptive RK45 integration of the geodesic equation for simplicity, spacetime flexibility, and handling of complicated and highly dynamic radiative media (i.e. the 'slow-light' paradigm). This does cause it to perform comparatively poorly to other ray tracers in situations with highly symetric spacetime and simple objects (e.g. a Schwarzschild or Kerr black hole with a thin accretion disk), so don't expect real time rendering.
* GPU integration for ~15x speedup (always looking for more). Currently only works with CUDA, but I intend to make it work on other platforms.
* Skyboxes/Textures can also be used on boundaries (e.g. with the wormhole).
* Synchrotron radiation code, resulting in the first examples I've seen of true color images rendered with a physically-based synchrotron spectrum. (I'm sure someone's done it, but I haven't seen one).
* Currently implemented scenes are: a Sun-like sphere in flat spacetime (spherical and cartesian) coordinates, Schwarzschild and Kerr black holes with constant temperature thin disks, and the wormhole from *Interstellar*.

## WIP:
* Ray-bundles / gaussian anti-aliasing: can use ray bundles like in *Interstellar* (James et. al. 2015, https://arxiv.org/abs/1502.03808) for anti-aliasing the skybox. My implementation currently ignores relativity, because the ray bundle equations are very complicated and probably extremely computationally expensive. Which ironically just makes it a rather ineficient anisotropic gaussian antialiasing scheme. Only used for the wormhole and the cartesian flat space test.
* A realistic model, based on state-of-the-art research, of a black hole along with it's accretion disk and jets that's also suitable for this ray tracing code. I chose to model M87*, also known as Pōwehi, due to it being the best observed black hole at the event horizon scale and known to have a jet. So far, I have crude model that is very far from what observations show (hint: not a thin-disk and cone), and kind of far from simulations as well (not a thin-disk). It does have synchrotron radiation and magnetic fields, though.
* The preprocessing and postprocessing of images is currently very inefficient, especially when using textures and ray bundles. Lots of room for improvement.

## Planned Features:
* A dedicated cameras.jl file that will include setups for rendering batches of images (e.g. for video). This could also help with the pre- and post- processing bottlenecks.
* Neutron star models; quiescent, pulsar, and accreting.
* Distributed computing capability (e.g. for video).
* Toggleable ability to output raw spectral brightness data (i.e. something like FITS files).
* A GUI
* Better documentation, including a nice PDF that is both a description of implementation and a manual.
* Ability to swap scenes without restarting Julia.

## Possible Future Features:
* Pre-compilation or maybe even a .exe for Windows.
* Hamiltonian-based numerical integration scheme, possibly activating automatically when it suits the situation.
* Kerr-Schild coordinates for a black hole.
* Taking rays from previous runs and continuing (e.g. rendering close things in approximately flat space, then switching to proper GR for everthing else).
* Making wacky stuff like geometrically closed universes.
* Black hole collisions from a lookup table metric.
* Rewriting in C++ or Fortran.

# Pretty Pictures
* A DNEG wormhole. Near side is the Milky Way, far side is artist's impression of a red nebula. ![Edited description from Gemma 4: A rendering of a wormhole dominates the center of the frame against a dark background of outer space. The sphere is characterized by swirling, turbulent patterns of deep crimson and bright red clouds, giving it a textured, gaseous appearance. A distinct, glowing red ring or halo borders its circumference, creating a sharp contrast with the surrounding space. The background is a vast expanse of black, densely peppered with tiny, white points of light representing distant stars. Wisps of dark, dust-like nebulae or cosmic clouds stretch diagonally across the frame behind the central sphere, adding depth to the composition. The overall mood is one of cosmic mystery and grand scale.](main/Assets/wormhole.png)
* A view from a r = 5M almost-circular orbit, above the disk of the dummy-model Kerr black hole, which has been thickened to be just barely optically thick for dramatic effect. ![Edited description from Gemma 4: A zoomed-in, macro view of a black hole simulation, showcasing the extreme distortion of light near an event horizon. The image features a warm color palette of creamy whites, soft tans, and golden-brown hues against a pitch-black background. A brilliant, sharp curve of blue-white light sweeps down from the top center and bends sharply toward the right. In the lower half, the light transforms into a softer, more diffused golden haze that stretches horizontally under the camera.](main/Assets/kerr_denser_disk.png)
* A view from a high-ish orbit of the Kerr black hole. ![Edited description from Gemma 4: A pitch-black background. A bright, hazy white accretion disk flows horizontally across the center, where it's image is warped by gravity to form a luminous, circular ring around a central dark sphere. The light is most intense and concentrated near the center, gradually fading into a soft, ethereal glow as it stretches outward to the left and right. The light has a blueish color on the left side, and a redish color on the right.](main/Assets/kerr_far.png)
* A wide-angle (135 degrees) view from inside the disk in a r = 3M orbit. ![Edited description from Gemma 4: A wide angle view of a black hole's gravitational lensing effect. A horizontal beam of light passes through a large, central black void, where it is dramatically warped into a bulbous, curved shape. The light glows with a cool, icy-white tone with a subtle bluish tint, contrasting sharply against the solid black of the space and the event horizon. The image is perfectly symmetrical, showing the light stretching and curving outward to form a luminous, ethereal halo around the core.](main/Assets/kerr_135deg_3M_in.png)
* A view from a circular orbit inside the ergosphere of the black hole, somewhat near the horizon. The camera has been rotated couterclockwise 90 degrees, so the accretion disk points up. ![Edited description from Gemma 4: A simulation of a black hole's event horizon from a perspective where the horizon sits at the bottom of the frame. The light features a cool, icy-blue and white color palette against a deep black background. A vertical pillar of light rises from the void, in the center of the fram, flaring outward as it goes. This is intersected by a wide horizontal band of light that curves dramatically upward at both ends, creating a symmetrical, bowed shape.](main/Assets/kerr_ergosphere.png)
* A view of my M87* model, near the jet. The jet is blue because of the synchrotron radiation. ![Edited description from Gemma 4: The image features a warm color palette of cream, beige, and soft brown. On the left, a large, hazy halo of light surrounds a central dark sphere. From the right side of this core, a brilliant blue-white, conical jet of light projects horizontally toward the right edge of the frame. The jet is most luminous at its origin point and gradually widens and fades into the surrounding pitch-black void.](main/Assets/simple_jet.png)

# Running it Yourself
### Requirements:
* Julia command line interface (e.g. standard REPL) and some packages (see Installation).
* It's only been tested in Ubuntu Linux. But it might work in Windows, since Julia runs in Windows.
* If you don't have a CUDA capable GPU, read the comments in Packages.jl to avoid trying to install CUDA.jl .
### Installation:
1. Install Julia and your REPL.
2. Clone the RTWCA Repository with git.
3. If you want to use something with textures, like the wormhole, you need to manually download the textures folder, which is seperate from the REPO to avoid bloat. Go to this Google Drive link https://drive.google.com/drive/folders/1AUyp-tCAWDCX97zUNzBwM4ZUSKzG_PZy?usp=sharing and move the "Textures" folder to sit inside the "main" folder.
4. Inside the "main" folder, activate Julia with whatever number of threads you want (on Linux with the Julia REPL, I run "julia -t auto" to automatically use all my CPU cores)
5. Run "include("Packages.jl")" to install all necessary packages (read the comments in that document about disabling CUDA and test mode if you want to remove those). Then wait like half an hour for it to compile the binaries for them all (don't worry, it only needs to do this once)
### Rendering Your First Image: a Kerr Black Hole
1. Run "include("initialization.jl")" and type "cpu" or "gpu" on the first prompt . On the second prompt type anything but "yes" (this prompt isn't case-sensitive, so you shouldn't type "YES" either).
2. Type:

        position, direction, pointing, β = initialize_world("kerr_boyer_lindquist")

to load the code for the Kerr scene and load the default camera parameters for this scene.
3. Choose how densly to sample the visible spectrum. I recomend 30~60 nm as a good balance between performance and accuracy:

        colors = range(350,step=30,stop=750)

4. Choose your horizontal number of pixels. I recommend 100~200 for the CPU, and 500~2000 if you have a compatible GPU:

        hor_pix = 500

3. To run, type: 

        img = gen_image(camera_pos = position, colors = colors, camera_dir = direction, max_dt_scale = 1e-1, max_steps = 1e4, x_pix = hor_pix, speed = β, camera_point = pointing);

This might take awhile to run on the CPU, and less time to run on the GPU, but with a long upfront compile time (~30s for me). Future runs should require little to no compilation.
It will return your image as an RGB array named "img" (or whatever you put there)
4. To display the image, type "imshow(img);" and type "save("image_name.png",img)" to save the image as a png.
5. To use a different scene, restart julia and give a different name to "initialize_world". To exit Julia, type "exit()".

# LLMs
I've started experimenting with using LLMs for some of the coding. So far, they've only been used for writing scripts and refactoring the Christoffel symbol functions. I don't think I'll ever get much past the point of having them write small pieces of code. Also, I used one to write most of the alt-text for the images in this README; Gemma 4 seems to write better sounding image descriptions than I do, and much faster. The rest of the README is all me.