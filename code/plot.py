from whole import *
from qutritalgo import *

x=[4, 6, 8, 10, 14, 18, 24, 32, 44, 58, 80, 106, 144, 194, 262, 354, 476, 642, 866, 1168, 1574, 2122, 2860, 3858, 5200]
y1=[0.3300023701419617, 0.4352737125216861, 0.4627548004106815, 0.35283483252731523, 0.3508257475900742, 0.3870241512512928, 0.30978864553710134, 0.2324395598261591, 0.17257997073841924, 0.24080280884214336, 0.28078571538128455, 0.2767748072900079, 0.22163581385876094, 0.2369755474126074, 0.2556590141691308, 0.17009362049627902, 0.2056224783718957, 0.18177678446227816, 0.1836447529697363, 0.14765976952611495, 0.11224231928106598, 0.12557592059204914, 0.144247085527759, 0.17721654353443334, 0.1854520638955331]
var1=[0.08820945642588411, 0.12740620573340858, 0.12520486375781192, 0.10713215978546814, 0.13513444415848191, 0.1114471198241485, 0.05276479368637146, 0.05894090609722423, 0.021521514722828355, 0.06570614543434489, 0.08038928675174062, 0.07603575831764119, 0.056925076749928506, 0.04939423982932273, 0.0468845281393811, 0.04939248741518238, 0.026922305652254393, 0.015524572510435207, 0.06889946745543614, 0.018206491062463023, 0.008855558286143188, 0.01792109259762013, 0.02211976533400147, 0.0658185313983073, 0.01150644709813845]
amax1=[0.9009433432609343, 1.2366762249679202, 0.9505863401384762, 1.271206637512612, 1.1114127974684456, 1.0827374213267873, 0.7776858591747041, 1.0769385935170628, 0.6342304167499211, 0.9355578916819163, 0.9726314941780465, 0.8402580021666456, 1.0426420031922081, 1.0122209881275044, 0.5991520287735467, 0.9302992476226599, 0.6466995124429549, 0.39978941259058487, 1.044320145612707, 0.5098280607256921, 0.3727010737018529, 0.46297757991544963, 0.3982510369513784, 1.1350262970508378, 0.33828783274350577]
y2=[0.22348966383349905, 0.26386738081543254, 0.2297011468265962, 0.32419815014624037, 0.23648254790260528, 0.17522427610228614, 0.13977080164298591, 0.2079458462463167, 0.23445170820041877, 0.14321606070305165, 0.18175456036731996, 0.09842773899762322, 0.18093451656978016, 0.1390209616115248, 0.1664739561505182, 0.15525720127030754, 0.15486077522700828, 0.17013239122921897, 0.08045355640179078, 0.08211327949093292, 0.07884799228191912, 0.07339275582394673, 0.09743304150573712, 0.0730016124948514, 0.08135436325761349]
var2=[0.025220404157455325, 0.03927471945343473, 0.02954940648571907, 0.05486450792658422, 0.06109308020940897, 0.01929504737192458, 0.020024761123907296, 0.026569078235754183, 0.03397781859811795, 0.011016864933657833, 0.030152049756366966, 0.008866723726860194, 0.03209465282367238, 0.02536628542472831, 0.04077336981944772, 0.03840648175439232, 0.027678907105903428, 0.038759293186369864, 0.0033975972073157818, 0.005929434187745195, 0.006433468205789153, 0.008901194793274983, 0.005529527439101883, 0.006720309204251536, 0.007923691831771136]
amax2=[0.5102655036144668, 0.7793372242699925, 0.6168871315774357, 0.8334665142479938, 0.8324144752481024, 0.4490645911141935, 0.5463894381228477, 0.5249452089992108, 0.5405891355116825, 0.3848434202345337, 0.6777831679406622, 0.4196171039469771, 0.6247261789000619, 0.6165710965043321, 0.7030801367324548, 0.8152340501616068, 0.5746550423897155, 0.7344169207995025, 0.19094470787683915, 0.2861933343003523, 0.2739243588249484, 0.42388660595275923, 0.2831543115699632, 0.32450753810249927, 0.3102211891008656]
y3=[0.3536184446374616, 0.5123666277536557, 0.5201574562361002, 0.31486874933499454, 0.3807519632957427, 0.2771316198187798, 0.27423370618747017, 0.27297244077103144, 0.2227075609772144, 0.3083901028013, 0.3202929365973971, 0.21723096680878987, 0.16658777195046254, 0.20897342578767536, 0.2194697027613509, 0.22995715221094054, 0.16107682575766108, 0.2806171924442324, 0.13374042521914553, 0.14202635144546047, 0.2110585169407171, 0.21343840939238098, 0.13700138580921858, 0.13450159530602562, 0.1260658522177982]
var3=[0.1356841522200855, 0.14652141107250993, 0.10927460016070029, 0.10031287801243664, 0.17122131876534052, 0.07133751000153082, 0.0899095467759677, 0.07414327910672462, 0.026751474714852366, 0.13382079694444177, 0.15852092089660155, 0.043454885228825874, 0.04964692530788133, 0.052074715894255795, 0.07842477323779624, 0.11579555830460714, 0.08590672347811157, 0.08845499592712046, 0.020753915443274495, 0.014569379744880633, 0.0353912661442819, 0.04378102960923398, 0.020309340876408503, 0.015497443121267064, 0.012878854172287709]
amax3=[1.1784880362708554, 1.1741238749928822, 1.10492055656272, 1.1404436762070935, 1.1925358996608104, 0.966075913137222, 0.8777725878894919, 1.201146330575572, 0.5120060400794647, 1.3343927871486982, 1.3300029927603862, 0.793362785569268, 0.8145165183430471, 0.9670741458483396, 1.0996519058594612, 1.3945087511040841, 1.3042768910275175, 1.1771032420015985, 0.5233148194981271, 0.4483058266722056, 0.6265842473952307, 0.6809978170625391, 0.4974364916107428, 0.3996527055139428, 0.3427760030847548]
y4=[0.21517248708861594, 0.16585208895500828, 0.1951951617985408, 0.2575238824187561, 0.16864282341254247, 0.11102755328406891, 0.22568763519132928, 0.12430188161654396, 0.13856369413194725, 0.17672017925447608, 0.20609541124162037, 0.16258785036098072, 0.11100969530459726, 0.11151935257965315, 0.07707119531826932, 0.13352987081325413, 0.06816766794938368, 0.07364133130321521, 0.05965041709180826, 0.06863333760904293, 0.09082681006467337, 0.06098249240151592, 0.0918094279987621, 0.05505354719117014, 0.08414529858454178]
var4=[0.038332324271514, 0.020376050158695383, 0.03911952122964173, 0.04282528522887342, 0.03141211336443539, 0.01791058034723929, 0.04741044285731936, 0.017323913444372674, 0.0140388227353651, 0.02755178373151583, 0.04172185180766346, 0.01742724637715987, 0.0077726703402669715, 0.025528346550530483, 0.007205767213754101, 0.024710743980228263, 0.00413236899999489, 0.008798169923673924, 0.002484685857271547, 0.0026588144609191724, 0.011657031459495032, 0.003236615127747685, 0.018448633921756687, 0.004133741440915083, 0.009858777777917497]
amax4=[0.610532065416516, 0.5584195589239277, 0.6273702406854664, 0.7099806307829779, 0.6471961551713767, 0.5966307152683327, 0.63382996039832, 0.46407718680234555, 0.4975899585146318, 0.4980950038572415, 0.6487904082223338, 0.45592305630227387, 0.2663629941215861, 0.6510753528838187, 0.30892826813495583, 0.6297511408839456, 0.21209161839618165, 0.309395789991744, 0.15486754204769468, 0.22002936569506917, 0.38083167093978193, 0.19553482203803266, 0.5326205198447784, 0.29471604139983515, 0.411559940765437]
fig, axs =plt.subplots(nrows=2,ncols=2,sharex=True)
axs[0, 0].errorbar(x,y1,yerr=var1,fmt='o',elinewidth=1,capsize=3,lw=0)
axs[0, 0].set_title("lambda = 1/2 ")
#axs[0, 0].scatter(x,amax1,c='r',marker=".")
axs[1, 0].errorbar(x,y2,yerr=var2,fmt='o',elinewidth=1,capsize=3,lw=0)
axs[1, 0].set_title("lambda = 2")
#axs[1, 0].scatter(x,amax2,c='r',marker=".")
axs[0, 1].errorbar(x,y3,yerr=var3,fmt='o',elinewidth=1,capsize=3,lw=0)
axs[0, 1].set_title("lambda = 1/4")
#axs[0, 1].scatter(x,amax3,c='r',marker=".")
axs[1, 1].errorbar(x,y4,yerr=var4,fmt='o',elinewidth=1,capsize=3,lw=0)
axs[1, 1].set_title("lambda = 4")
#axs[1, 1].scatter(x,amax4,c='r',marker=".")
plt.xlabel('n')
#plt.ylabel('ratio')
plt.xscale('log')
#plt.title('Transverse field Ising model')
plt.show()