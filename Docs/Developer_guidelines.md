# SIMsalabim project: Developer Instructions

All input that can help us make SIMsalabim better and more useful for users is welcome. This document contains instructions for what to do if you would like to contribute to or tinker around with the code.

## Filing Issues
If anything is unclear or you have any feature requests, you can file an 'Issue' in the Issues section on GitHub. Please don't hesitate to do this, as it might be of help for other user or us as maintainers as well. If you have a question we can provide help in the issues section so that future users can use if for reference. 

If you have any feature requests, we would first need to discuss the how and why in detail, which we would also do in the issues section. We will decide on a plan of action on a per case basis, but expect to make many of the code changes ourselves. 

In case of bug reports we will need the exact configuration of the 'device_parameters.txt' file as many appendant bugs can be traced back to settings that result in an unphysical device. Including that, along with the terminal output (including version of the program used), will help us dial in on the problem quickly.


## Overall structure
- ZimT and SimSS share much of the same code, so try to put new code in the units, not in ZimT or SimSS
- Routines, types and constants that are specific to drift-diffusion modelling should go in DDRoutines or DDTypesAndConstants. These units should have the same version number as the ZimT/SimSS codes.
- Somewhat generic routines, types and constants should go into units TypesAndConstants, InputOutputUtils, and NumericalUtils. As yet, these units do not have a version number.
- Exit Codes are defined in unit TypesAndConstants.
- Input parameters and parameters that are directly derived from those input parameters (like Booleans based on a 0,1 input) are stored in par (='parameters') of type TInputParameters. The order of the fields within this record (integers, reals, etc are grouped) should follow the order of the variables in the file device_parameters.txt. Note: some parameters are limited to either SimSS or ZimT. However, we do not make use of the possibility of introducing a variant part in these records as this led to unexpected behaviour. For example, when a field is only defined (through the variant part) to ZimT, it was still availalbe in SimSS and no warning or exception was given when using that part of the record (with an undefined value!).
- Other variables that remain do not change throughout the simulation are stored in stv (='static variables') of type TStaticVars.

## General remarks on changing the code
- Make sure you understand what is permitted&mdash;and what is not&mdash;under the GLPL licence.
- Describe changes in file change_log.txt: Briefly describe the changes to the code, but this is not the main point (as one can easily figure them out using a difference viewer). More importantly, motivate and/or explain why the change was made.
- Run the tests as outlined in the folder 'Tests'. These are just a few tests to assess the basic functionality so additional testing is required.
- Magic numbers used in the code are, mostly, listed as constants in DDTypesAndConstants. Here we also indicate where&mdash;in which routine&mdash;the magic number is used and what it does.
-Try to limit local variables in and parameters passed to subroutines to 32K as this is the limit of some processors. 

## Notation within the codes
- Use Linux line-endings
- Add ample comments to the code, using curly brackets {}. Describe *why* rather than *what* when formulating comments.
- Naming of variables that represent a physics unit: snake-case where multiple sub- and superscripts are separated by an underscore, for example: E_t_bulk
- Naming of other variables: lower camel-case, i.e. camel-case where the first word is always lowercase, for example: bulkTrapType
- Naming of procedures and function: pascal case combined with snake-case: Calc_Elec_Mob, for example
- FPC/Pascal key words fully capitalized: IF ... THEN ... ELSE
- Indentation width: 4

## Parameters
- Physical units: SI only, except for work functions (eV). Work functions are given as the distance to vacuum so they are positive.
- Device parameters: try to keep the same order for ZimT and SimSS. Every parameter should have a unit (if any) and a short description. After reading the parameters, they should be checked (to some extent) whether they make sense (logical, physical) in Check_Parameters.
