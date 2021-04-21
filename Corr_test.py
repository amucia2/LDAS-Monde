

r_decid_c = correlation(lai_cgls,vodc,patch_3,perc)
r_decid_win_c = correlation(lai_cgls_win,vodc_win,patch_3,perc)
r_decid_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_3,perc)
r_decid_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_3,perc)
r_decid_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_3,perc)

r_conif_c = correlation(lai_cgls,vodc,patch_4,perc)
r_conif_win_c = correlation(lai_cgls_win,vodc_win,patch_4,perc)
r_conif_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_4,perc)
r_conif_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_4,perc)
r_conif_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_4,perc)

r_c3_c = correlation(lai_cgls,vodc,patch_6,perc)
r_c3_win_c = correlation(lai_cgls_win,vodc_win,patch_6,perc)
r_c3_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_6,perc)
r_c3_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_6,perc)
r_c3_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_6,perc)

r_grass_c = correlation(lai_cgls,vodc,patch_9,perc)
r_grass_win_c = correlation(lai_cgls_win,vodc_win,patch_9,perc)
r_grass_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_9,perc)
r_grass_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_9,perc)
r_grass_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_9,perc)

perc = 0.1
r_veg_c = correlation(lai_cgls,vodc,patch_nonveg,perc)
r_veg_win_c = correlation(lai_cgls_win,vodc_win,patch_nonveg,perc)
r_veg_spr_c = correlation(lai_cgls_spr,vodc_spr,patch_nonveg,perc)
r_veg_sum_c = correlation(lai_cgls_sum,vodc_sum,patch_nonveg,perc)
r_veg_aut_c = correlation(lai_cgls_aut,vodc_aut,patch_nonveg,perc)


headers = ['Veg Type','All Seasons','Winter','Spring','Summer','Autumn']
decid_c = ['Deciduous',r_decid_c[0],r_decid_win_c[0],r_decid_spr_c[0],r_decid_sum_c[0],r_decid_aut_c[0]] 
conif_c = ['Coniferous',r_conif_c[0],r_conif_win_c[0],r_conif_spr_c[0],r_conif_sum_c[0],r_conif_aut_c[0]] 
c3_c = ['C3 Crops',r_c3_c[0],r_c3_win_c[0],r_c3_spr_c[0],r_c3_sum_c[0],r_c3_aut_c[0]]
grass_c = ['Grasslands',r_grass_c[0],r_grass_win_c[0],r_grass_spr_c[0],r_grass_sum_c[0],r_grass_aut_c[0]]
datac = [decid_c,conif_c,c3_c,grass_c]
cortable_c = pd.DataFrame(data=datac,columns=headers)