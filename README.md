# BPEM
A Barotropic Primitive Equation Model
This model was written for study the atmosphere model, it might be a tool for researching the barotropic atmosphere.

#Install the model
Adjust the configure file configure.BPEM to fitting your environment, then type
make
If the compile finish sucessful, the executable BPEM.exe will be shown in the root directory of BPEM

#Input source data and initialization

This model needs a ERA-interim netCDF file to provide the first guess field(FG), and the lateral boundary condition(LBC).
You may download ERA-interim data from "http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=pl/"
Required fields are "Geopotential","U component of wind" and "V component of wind" at 500hPa(for original 
barotropic model).
After downloading the ERA-interim netCDF file, you may generate the FG and LBC by using the Matlab script "initialize/real_em.m".
In "real_em.m", I support 2 methods to generate the model domain
1. If you have a WRF Preprocess System(WPS) geogrid output file "geo_em.d01.nc", then you may generate the domain by using it,
of course the domain will be the same as which in WRF.
2. Also, you can generation the domain by input the lambert projection parameters, you might see the parameter in real_em.m

Here I've put some test data in BPEM/input
geo_em.d01.nc is the WPS geogrid output.
201710_ZUV.nc is the ERA-interim file.

#Namelist

Once you've generated the FG and LBC, for letting the model to run as the way you want, you should edit the namelist.input,
system_type         : choose from 'Windows' or 'Linux'
max_openmp_threads  : the maximum openmp threads number you want to use, it depends on how many threads do your computer has.
run_hours           : how many hours do you want to run the model, if you're using the test data I provided, 720 is the max value.
history_interval    : how ofter do you want to output the forecasting data, unit is minute.history_interval*60 should be a 
                      integer multiple to time_step, or the model won't output correctly.
input_data_path     : where the input files are, for the test case ./BPEM/initialize
output_data_path    : where do you want to store the output data, for the test case ./BPEM/output
time_step           : the integration time step, if the time step is too large, the model might blow up. for 100km resolution,
                      200s is a good choice.
integration_option  : Choose from 1 for loop integration, or 2 for matrix integration, but for now matrix scheme is not avaliable.
z0                  : not used, please set to 0
smooth_coefficient  : not used
spec_exp            : exponential multiplier for the relaxation zone ramp, used with a specified boundary condition. 0. = linear ramp,
                      default; 0.33 = ~3*dx exp decay factor. This may be useful for long simulations.

#Postprocess

BPEM will output a netCDF file, you may plot the data by some useful software, such as "ncview","Panoply","NCL","Matlab", here I've
wrote a simply Matlab script in ./BPEM/postprocess/post_process.m.

If you have any question about this model, please send e-mail to : dwyanechou@qq.com
