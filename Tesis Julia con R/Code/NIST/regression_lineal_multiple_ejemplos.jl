using DataFrames, GLM, CSV

# Primer ejemplo
data_kid = CSV.read("C:/Users/Valeria/Documents/ITAM/Tesis/Julia con R/Regression_and_other_stories/ROS-Examples-master/ROS-Examples-master/KidIQ/data/kidiq.csv", DataFrame)

fm = @formula(kid_score ~ mom_hs + mom_iq + mom_hs*mom_iq)

kidscore_lm = lm(fm, data_kid)

# Segundo ejemplo
# No funciona por los NA
earnings_data = CSV.read("C:/Users/Valeria/Documents/ITAM/Tesis/Julia con R/Regression_and_other_stories/ROS-Examples-master/ROS-Examples-master/Earnings/data/earnings.csv", DataFrame)
earnings_data.cHeight = earnings_data.height .- 66

fm = @formula(weight ~ height + male + ethnicity)
earnings_lm = lm(fm, earnings_data)

#Segundo intento con CSVFiles
earnings_data = DataFrame(load("C:/Users/Valeria/Documents/ITAM/Tesis/Julia con R/Regression_and_other_stories/ROS-Examples-master/ROS-Examples-master/Earnings/data/earnings.csv"))
earnings_data.cHeight = earnings_data.height .- 66

fm = @formula(weight ~ cHeight + male + ethnicity)
earnings_lm = lm(fm, earnings_data)
