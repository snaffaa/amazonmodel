"""
amazon_sediment_production_functions: collection of functions - all taking direct input and 
returning output- that are intended to serve as componenents(monthly soul loss)of the 
Amazon sediment Transport Model (PCR-GLOBWB-erosion and sedimentation module_.

"""

# modules
import pcraster as pcr
import os
outdname ='/scratch/safaa/output/1980/new_cover_fraction_run'

def compute_slope_length_factor(slope_length, slope_gradient):
    '''
compute_slope_length_factor: function that returns the factor that corrects the 
erosion for the slope length. The function uses the relationship used by Renard 
et al. (1997) where the power gamma is fitted to a continuous function, capped 
unity (when the slope angle is 32 degrees or the slope 0.625).

    Input:
    ======
    slope_length:           length of the contributed slope projected on the
                            horizontal, in [m];
    slope_gradient:         gradient of the contributing slope, in [m/m].

    Output:
    =======
    slope_length_factor:    correction factor to take the slope length into
                            account [1].

'''

    # compute the power gamma first, then return the slope factor
    gamma = pcr.min(1.00, 1.183156498 * slope_gradient ** 0.351920534)
    
    slope_length_factor = (slope_length / 22.13) ** gamma
    
    # return the factor
    return slope_length_factor
    
#fed

def compute_slope_steepness_factor(slope_gradient):
    '''
compute_slope_steepness_factor: function that returns the factor that corrects 
the erosion for the slope steepness. The function uses the relationship used by 
Nearing (1997).

    Input:
    ======
    slope_gradient:         gradient of the contributing slope, in [m/m].

    Output:
    =======
    slope_steepness_factor: correction factor to take the slope steepness
                            into account [1].

'''

    # compute the slope angle from the gradient; this uses the default option
    # of PCRaster (e.g., radians)
    slope_angle = pcr.atan(slope_gradient)
    slope_gradient= slope_gradient *100

    slope_steepness_factor = -1.5 + (17.0 / (1 + \
            pcr.exp(2.3 - (6.1 * pcr.sin(slope_angle)))))
    
    # return the factor
    return slope_steepness_factor
    
#fed

def compute_slope_factor(slope_length, slope_gradient):
    '''
compute_slope_factor: function that returns the factor that corrects the 
erosion for the combined effect of slope length and slope gradient.
The function uses the relationship used by Renard 
et al. (1997) where the power gamma is fitted to a continuous function, capped 
unity (when the slope angle is 32 degrees or the slope 0.625), multiplied with 
the function for the slope steepness developed by Nearing (1997).

    Input:
    ======
    slope_length:           length of the contributed slope projected on the
                            horizontal, in [m];
    slope_gradient:         gradient of the contributing slope, in [m/m].

    Output:
    =======
    slope_factor:           correction factor to take the slope length and the 
                            slope steepness into account [1].

'''

    # return the factor
    return compute_slope_length_factor(slope_length, slope_gradient) * \
            compute_slope_steepness_factor(slope_gradient)
    
#fed

def compute_conservation_factor(cover_fraction_info, P_factor_info, selected_slope):
    '''
compute_conservation_factor: The function uses the global estimation of Pham
et al.(2003)for every land cover type expect crop lands. For crop lands the
function calculate the conservation factor on the basis of the slope gradient
using the global method of Wener (1981).

    Input:
    ======

    selected_slope:         gradient of the corrosponding slope of 
                            crop land types in [m/m].
    
    conservation_factor_info:land cover types and teir estimated 
                            conservation factor values
    Output:
    =======
    conservation_factor:    factor [1] that corrects for the presence of soil
                            conservation measures (1 if none).) 

'''
    P_factor = pcr.scalar(0)
    # iterate over the cover fractions
    for key in cover_fraction_info.keys():
        if key == "croplands":
            P_factor_covertype =  0.2 + 0.3 * selected_slope 
        else:
            P_factor_covertype =   P_factor_info[key]
       #fi
        P_factor = P_factor \
                + cover_fraction_info[key] * P_factor_covertype

    #rof
    
    return P_factor 
#fed

def compute_interception_fraction(cover_fraction_info, interception_info):
    part_int0 = pcr.scalar(0)
    # iterate over the cover fractions
    for key in cover_fraction_info.keys():
        part_int0 = part_int0 \
                + cover_fraction_info[key] * interception_info[key]
    #rof
    part_int1 = pcr.max (0 , part_int0)
    part_interception = pcr.min (1, part_int1)
    
    return part_interception
#fed

def compute_annual_erosivity(rain_annual):
    '''
compute_annual_erosivity: function to compute the annual erosvity based on the 
relationship of Renard and Freimund (1994).

    Input:
    ======
    rain_annual:        annual rain in [mm].
    
    Output:
    =======
    eros_annual:        annual erosivity, R-factor [MJ * mm / ha /hr / yr].

'''

    # compute the erosivity
    return pcr.ifthenelse(rain_annual <= 850.0, \
            0.0483 * rain_annual ** 1.61, \
            587.8 - 1.219 * rain_annual + 0.00415 * rain_annual ** 2)

#fed

def partition_erosivity(rain_month, rain_annual, eros_annual):
    '''
partition_erosivity: function to partition the annual erosvity into smaller 
parts based on the rainfall.Effective rainfallreplaced the monthly
rainfall and calculated according to Morgan (2005)


    Input:
    ======
    rain_month,
    rain_annual:        part and total annual rain in [mm].
    part_interception:  proportion of rainfall interception by the vegetation 
                        or crop cover between 0 and 1
    eros_annual:        annual erosivity, R-factor [MJ * mm / ha /hr / yr].

    
    Output:
    =======
    eros_part:          part of the erosivity.

'''
    # compute the erosivity
    return pcr.ifthenelse(rain_annual > 0, \
        rain_month / pcr.max(rain_month, rain_annual), 0) * eros_annual
#fed

'''
compute_vegcover_factor: function to compute the erosion reduction factor due 
to the vegetation cover (crop management factor) using vegetation fraction 
and annual erosivity by Xihua Yang (2014)
    Input:
    ======
    ground_cover_fraction:  GCj monthly ground cover calculated by adding the maps 
                            of cover fraction from all land cover types.                
    
    
    annual ersovity index: AEIj calcaluated from dividing R(monthly Erosivity),
                                    to R (annual erosivity)
                            AEj = R_month/R_year
    
    Output:
    =======
    vegcover_factor:        vegetation cover (crop management)factor [1].
'''
def compute_vegcover_factor2(ground_cover_fraction,R_month,R_year ):
    AEI = R_month/ R_year
    gc_tmp = -0.799 - 7.74 * ground_cover_fraction
    pcr.report(gc_tmp, os.path.join(outdname, 'gc_tmp'))
    gcf_tmp=0.0449 * (ground_cover_fraction**2)
    pcr.report(gcf_tmp, os.path.join(outdname, 'gcf_tmp'))
    c_tmp= (gc_tmp + gcf_tmp)
    pcr.report(c_tmp, os.path.join(outdname, 'c_tmp'))
    c_ex_tmp=pcr.exp(c_tmp)
    c_ex_tmp= pcr.cover(c_ex_tmp, mv)
    pcr.report(c_ex_tmp, os.path.join(outdname, 'c_ex_tmp'))
    c_factor = c_ex_tmp * AEI
    c_factor= pcr.cover(c_factor, mv)
    pcr.report(c_factor, os.path.join(outdname, 'c_factor'))
    return c_factor

'''
compute_melt_time_fraction: function that returns the fractional time that the 
air temperature during a period (e.g., month) is below the the melt rate.

    Input:
    ======
    tmin:                 minimum air temperature (degC) in a period,
                          e.g., per month;
    tmax:                 maximum air temperature (degC) in a period;
    tmelt:                melt temperature (degC), default is 0 deg C.
    
    Output:
    =======
    melt_time_fraction:   time the air temperature is below the melt temper-
                          ature; the fraction is set to zero if the temperature
                          does not reach above freezing, it is set to unity if
                          the temperature is above freezing throughout.

'''
def compute_melt_time_fraction(tmin, tmax, tmelt = 0.0):
    # compute the fractional melt time
    melt_time_fraction = pcr.min(1.0, \
            pcr.max(0, tmax - tmelt) / pcr.max(0.1, tmax - tmin))

    return melt_time_fraction    
#fed

def compute_kfact_monthly(kfact_annual, kfact_seasonal, melt_time_fraction):
    '''
compute_kfact_monthly: function that computes the corrected, monthly erodibility 
factor based on the equation of Van Dijk (2001), including the time the 
temperature is above the melt temperature, the seasonality ratio that is texture 
based and the soil erodibility that applies on average over the year.
'''
    return (melt_time_fraction * (kfact_seasonal - 1.0) + 1.0) * kfact_annual
#fed

def return_default_soil_erodibility( \
        mass_fraction_sand, mass_fraction_silt, mass_fraction_clay,\
        organic_matter_content,\
        ):
    '''
return_default_soil_erodibility: function that returns the default soil 
erodibility as the function of the specified inputs using the modified soil 
erodibility computation of the USLE as presented by Torri et al. and with the 
logarithm of the mean particle size, Dg, replaced by an estimate based on the 
the textural fractions as presented by Yang et al. (2003).

    Input:
    ======
    mass_fraction_sand,
    mass_fraction_silt,
    mass_fraction_clay:     fractions by mass of the soil textural classes;
    organic_matter_content: organic matter content as percentage.
    
    Output:
    =======
    soil_erodibility:       the default soil erodibility
                            [tonnes * ha * h / ha / MJ / mm].

'''

    # median particle size after Yang et al. (2003)
    dg = -3.5 * mass_fraction_sand - 2.0 * mass_fraction_silt - \
            0.5 * mass_fraction_clay
    # organic matter divided by clay content
    omc = pcr.ifthenelse(mass_fraction_clay > 0.001, organic_matter_content / \
            pcr.max(0.001, mass_fraction_clay), 0.0)
    # get the exponential term
    e_omc = (-0.0021 * omc - 0.00037 * omc ** 2 - 4.02 * mass_fraction_clay + \
            1.72 * mass_fraction_clay ** 2)
    
    # computing the soil_erodibility
    soil_erodibility = 0.0293 * (0.65 - dg + 0.24 * dg ** 2) * pcr.exp(e_omc)
    
    # return the soil erodibility
    return soil_erodibility
#fed

def return_soil_erodibility_seasonality( \
        mass_fraction_sand, mass_fraction_silt, mass_fraction_clay):
    '''
return_soil_erodibility_seasonality: function that returns the ratio of the 
k-factor representing the soil erodibility that is linked to the seasonality;
values are based on the work by Hoch (2014).

    Input:
    ======
    mass_fraction_sand,
    mass_fraction_silt,
    mass_fraction_clay:     fractions by mass of the soil textural classes;
    
    Output:
    =======
    soil_erod_seasonality:  the default soil erodibility
                            [tonnes * ha * h / ha / MJ / mm].

'''

    # setting the seasonality ratio of the soil erodibility factor
    soil_erod_seasonality = pcr.scalar(1.44)

    # subdivide the soil into coarse, medium and fine on the basis of texture
    soil_erod_seasonality = pcr.ifthenelse(mass_fraction_clay > 0.35,\
            pcr.scalar(1.17), soil_erod_seasonality)
    soil_erod_seasonality = pcr.ifthenelse((mass_fraction_clay < 0.18) & \
                (mass_fraction_sand > 0.65),\
            pcr.scalar(4.50), soil_erod_seasonality)
    soil_erod_seasonality = pcr.ifthenelse((mass_fraction_clay <= 0.35) & \
                (mass_fraction_sand <= 0.65),\
            pcr.scalar(1.44), soil_erod_seasonality)  

    # return the seasonality of the soil erodibility
    return soil_erod_seasonality

#fed

def compute_soil_loss( \
        erosivity_period, erodibility_period, \
        slope_factor, vegcover_factor, conservation_factor, CellArea):
    '''
compute_soil_loss: function that computes the soil loss according to the RUSLE 
method of Renard et al. (1997) as the product of the input factors.

    Input:
    ======
    erosivity_period:       erosivity per period, R-factor in units of
                            [MJ * mm / ha /hr / yr];
    erodibility_period:     the default soil erodibility
                            [tonnes * ha * h / ha / MJ / mm] for the period;
    slope_factor:           correction factor to take the slope length and the 
                            slope steepness into account [1];
    vegcover_factor:        vegetation cover (crop management) factor [1].
    conservation_factor:    factor [1] that corrects for the presence of soil
                            conservation measures (1 if none).    

    Output:
    =======
    soil_loss:              soil loss according to the RUSLE method as the
                            product of the respective factors in units of
                            [tonnes per hectare per period].

'''
    # return the total soil loss in ton/year
    return erosivity_period * erodibility_period * \
            slope_factor * vegcover_factor * conservation_factor * \
            CellArea / 10000 
#fed

def compute_delivery_ratio(HydroCoeff, slope_gradient, Manning, slope_length):
    '''
compute_sed_transport: function that computes the sediment transport based on 
van Dijk work(2001) as the product of the input factors.

    Input:
    ======
    hydrological coefecient:   is relative value, proportional to runoff per 
                               period and rain per period [ mm / mm];
    slope_gradient:            the slope angle in percentage[%],
                             
    Manning's roughness coeffecient: explained in the compute_Manning function
                               [s / m^1/3];
    slope_length:        


    Output:
    =======
    delivery ratio:         delivery ratio according as the
                            product of the respective factors in units of
                            [].

'''   
    alpha = pcr.scalar(9.53)
    beta = pcr.scalar(0.79)
    slope_gradient_Perc = slope_gradient * 100

    tmp0 = (HydroCoeff * pcr.sqrt(slope_gradient_Perc)) / (Manning*slope_length)
    tmp1 = alpha *pow(tmp0 ,beta)
    tmp2 = pcr.max(0, tmp1)
    DeliveryRatio = pcr.min(1, tmp2)
    return DeliveryRatio
#fed

def compute_hydro_coeff(Qsurface, rain_month):
    tmp1 = Qsurface / rain_month
    tmp2 = pcr.cover(tmp1, 0.0001)
    tmp3 = pcr.max(0, tmp2)
    HydroCoeff = pcr.min(1, tmp3)
    return HydroCoeff
#fed


def Manning_n1(slope_gradient, Manning_n1_TBL):
    """
Manning_n1: function to compute a scalar PCRaster field with the 
Manning's n as a weighed average of slope angle percentage per cell and 
the corresponding value due to surface irrigularities.

Input:
======
slope_gradient:        Pcraster maps of slope fraction per cell;

Output:
=======
Manning_n1: 

Manning_v3: function to compute a scalar PCRaster field with the 
Manning's n as a weighed average as  result of changes in vegetation
cover per cell based on ndvi monthly maps 


Input:
======
ndvi:       normalized difference vegetation index. this index varies between 
            -1.0 and 1.0.

Output:
=======
Manning_v3:              
                       
"""
    slope_gradient_Perc = slope_gradient * 100
    return pcr.lookupscalar(Manning_n1_TBL, slope_gradient_Perc)
#fed

def Manning_v3(ndvi):
    return pcr.ifthenelse(ndvi < 0, 0, ndvi)
#fed


def Manning_n4(cover_fraction_info, manning_n_info):
    """
get_manning_n: function to compute a scalar PCRaster field with the 
Manning's n as a weighed average of the vegetation fraction per cell and 
the corresponding value.

Input:
======
cover_fraction_info:  dictionary with the land cover class as key 
                      and the vegetation fraction for that particular
                      vegetation type as a spatial scalar PCRaster
                      field;
manning_n_info:       dictionary holding the same keys as the 
                      cover_fraction_info and a float or a scalar
                      PCRaster field of the Manning's n.

Output:
=======
manning_n:              spatial scalar PCRaster field with the
                        cell-average Manning's n.
"""
    # set manning's n
    manning_n = pcr.scalar(0)

    # iterate over the cover fractions
    for key in cover_fraction_info.keys():
        # add the manning's n
        manning_n = manning_n \
                + cover_fraction_info[key] * manning_n_info[key]
    #rof

    return manning_n
#fed
    
def compute_Manning(Mn1, Mn4, ndvi):
    Mn2 = pcr.scalar(0.049)
    Mv3 = Manning_v3(ndvi)

    return Mn1 + Mn2 + Mv3 * Mn4
#fed
    
'''
compute_erosion: function that computes the sediment input to river based 
on Hoch (2014) as the product of the input factors.

    Input:
    ======
    soil loss:              [tonnes per period];      
                            
    delivery ratio:         [] for the period;

    Output:
    =======
    Erosion:               the amounts of sediment input to river 
                           [tonnes per period].

'''
def compute_erosion(soil_loss, DeliveryRatio):
    Erosion = soil_loss * DeliveryRatio    
    return  Erosion
#fed
