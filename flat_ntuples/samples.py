era_dict = {}

genericPath = '/gwteras/cms/store/user/fiorendi/p5prime/'

LMNR_sub_B  = 'DoubleMuonLowMass/B0KstarMM_Sept14_SoftMuHitsSkim_LMNR_ReReco23Sept_2016B/180914_100500/'
LMNR_sub_C  = 'DoubleMuonLowMass/B0KstarMM_Sept14_SoftMuHitsSkim_LMNR_ReReco23Sept_2016C/180914_100702/'
LMNR_sub_D  = 'DoubleMuonLowMass/B0KstarMM_Sept14_SoftMuHitsSkim_LMNR_ReReco23Sept_2016D/180914_100909/'
LMNR_sub_E  = 'DoubleMuonLowMass/B0KstarMM_Sept14_SoftMuHitsSkim_LMNR_ReReco23Sept_2016E/180914_101114/'
LMNR_sub_F  = 'DoubleMuonLowMass/B0KstarMM_Sept14_SoftMuHitsSkim_LMNR_ReReco23Sept_2016F/180914_101314/'
LMNR_sub_G  = 'DoubleMuonLowMass/B0KstarMM_Sept14_SoftMuHitsSkim_LMNR_ReReco23Sept_2016G/180914_101516/'
LMNR_sub_H  = 'DoubleMuonLowMass/B0KstarMM_Sept14_SoftMuHitsSkim_LMNR_ReReco23Sept_2016Hv2/180914_101757/'
LMNR_sub_H2 = 'DoubleMuonLowMass/B0KstarMM_Sept14_SoftMuHitsSkim_LMNR_ReReco23Sept_2016Hv3/180914_101959/'

##LMNR
# b: 0-2
# d: 0-1
# g: 0-1
# h: 0-1

PSI_sub_B  = 'Charmonium/B0KstarMM_Sept14_SoftMuHitsSkim_Charmonium_ReReco23Sept_2016B/180914_113156/'
PSI_sub_C  = 'Charmonium/B0KstarMM_Sept14_SoftMuHitsSkim_Charmonium_ReReco23Sept_2016C/180914_113400/'
PSI_sub_D  = 'Charmonium/B0KstarMM_Sept14_SoftMuHitsSkim_Charmonium_ReReco23Sept_2016D/180914_113959/'
PSI_sub_E  = 'Charmonium/B0KstarMM_Sept14_SoftMuHitsSkim_Charmonium_ReReco23Sept_2016E/180914_114201/'
PSI_sub_F  = 'Charmonium/B0KstarMM_Sept14_SoftMuHitsSkim_Charmonium_ReReco23Sept_2016F/180914_113601/'
PSI_sub_G  = 'Charmonium/B0KstarMM_Sept14_SoftMuHitsSkim_Charmonium_ReReco23Sept_2016G/180914_113802/'
PSI_sub_H  = 'Charmonium/B0KstarMM_Sept14_SoftMuHitsSkim_Charmonium_ReReco23Sept_2016Hv2/180914_112705/'
PSI_sub_H2 = 'Charmonium/B0KstarMM_Sept14_SoftMuHitsSkim_Charmonium_ReReco23Sept_2016Hv3/180914_112905/'
    
list2 = [
           '0000/B0ToKstMuMu_*.root',
           '0001/B0ToKstMuMu_*.root',
           '0002/B0ToKstMuMu_*.root',
]
list1 = [
           '0000/B0ToKstMuMu_*.root',
           '0001/B0ToKstMuMu_*.root',
]
list0 = [
           '0000/B0ToKstMuMu_*.root',
]



era_dict['LMNR_2016B' ] = [ genericPath + LMNR_sub_B , list2]
era_dict['LMNR_2016C' ] = [ genericPath + LMNR_sub_C , list0]
era_dict['LMNR_2016D' ] = [ genericPath + LMNR_sub_D , list1]
era_dict['LMNR_2016E' ] = [ genericPath + LMNR_sub_E , list1]
era_dict['LMNR_2016F' ] = [ genericPath + LMNR_sub_F , list0]
era_dict['LMNR_2016G' ] = [ genericPath + LMNR_sub_G , list2]
era_dict['LMNR_2016H' ] = [ genericPath + LMNR_sub_H , list2]
era_dict['LMNR_2016H2'] = [ genericPath + LMNR_sub_H2, list0]

era_dict['PSI_2016B' ] = [ genericPath + PSI_sub_B , list2]
era_dict['PSI_2016C' ] = [ genericPath + PSI_sub_C , list0]
era_dict['PSI_2016D' ] = [ genericPath + PSI_sub_D , list1]
era_dict['PSI_2016E' ] = [ genericPath + PSI_sub_E , list1]
era_dict['PSI_2016F' ] = [ genericPath + PSI_sub_F , list0]
era_dict['PSI_2016G' ] = [ genericPath + PSI_sub_G , list2]
era_dict['PSI_2016H' ] = [ genericPath + PSI_sub_H , list2]
era_dict['PSI_2016H2'] = [ genericPath + PSI_sub_H2, list0]

