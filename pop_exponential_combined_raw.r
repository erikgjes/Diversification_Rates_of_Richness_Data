

#pdf(file='./pop_exponential/literate_mcmc_logs/COMBINED_RTT_plots.pdf',width=12, height=8)
#par(mfrow=c(2,4))
#library(scales)
###COMBINED_mcmc###

#Number Shifts
sp_rate_pop_exponential_unique=c(2.0, 3.0,4.0,5.0,6.0)
sp_rate_pop_exponential_counts=c(237154, 154718,7831,285,12)
#sp_rate_ Plot
sp_rate_pop_exponential_time=c(-48.994897959183675, -47.984693877551024,-46.974489795918366,-45.964285714285715,-44.954081632653065,-43.94387755102041,-42.933673469387756,-41.923469387755105,-40.913265306122454,-39.9030612244898,-38.892857142857146,-37.88265306122449,-36.87244897959184,-35.86224489795919,-34.852040816326536,-33.84183673469388,-32.83163265306122,-31.821428571428573,-30.81122448979592,-29.801020408163268,-28.790816326530614,-27.78061224489796,-26.77040816326531,-25.760204081632654,-24.75,-23.73979591836735,-22.729591836734695,-21.71938775510204,-20.70918367346939,-19.698979591836736,-18.688775510204085,-17.67857142857143,-16.66836734693878,-15.658163265306122,-14.647959183673471,-13.63775510204082,-12.627551020408163,-11.617346938775512,-10.607142857142861,-9.59693877551021,-8.586734693877553,-7.576530612244902,-6.566326530612251,-5.556122448979593,-4.5459183673469425,-3.535714285714292,-2.525510204081634,-1.5153061224489832,-0.5051020408163254)
sp_rate_pop_exponential_rate=c(0.9807354711237168, 0.14195355732195283,0.14194847060092217,0.1419686024777027,0.1419930597980215,0.1420507054828937,0.1421084178487114,0.1421671756381956,0.14233857960325513,0.14239173939057456,0.14245585155193424,0.1425031119072479,0.1425818462415709,0.14267877490261277,0.1428043818889898,0.1428839011779715,0.1429477381055511,0.1430104739832358,0.14306596888914375,0.14314330980667817,0.14325008381443763,0.143409917372444,0.1437134368721448,0.14400923854733744,0.14423047181772108,0.14447538420239062,0.14471575352066268,0.14495861657457396,0.14527530538438999,0.14563748145447736,0.14614206363914814,0.14673733480998088,0.14746487668181174,0.14814371819392003,0.14899183467813437,0.14983213169444048,0.15078516820512286,0.15201691332428846,0.15374575681507693,0.15509804432431812,0.15595680710436843,0.15667727010360186,0.15718251544063436,0.15771280490645667,0.1584609800179472,0.15935286988053574,0.15965775161384754,0.15997909012203776,0.15998185346010882)
sp_rate_pop_exponential_minHPD=c(0.13024615813483595, 0.11310487507996313,0.11310803012572121,0.11346073101993867,0.11367563268784014,0.11382600940544897,0.11405945937909426,0.1143920701351628,0.11472427168112394,0.11450422526711065,0.11461497353973565,0.11473135663637288,0.11484617845721629,0.11517309963402046,0.11571226735246715,0.11580712044971885,0.11567731735134522,0.11581349910282192,0.11580649065210327,0.1158496443580042,0.11581375182713242,0.11624043007001673,0.11657405578471207,0.1167609641117092,0.11708074228112932,0.11709962439916263,0.11718326465253909,0.11720834846578747,0.1176500231705734,0.11764515443292371,0.11764502722520465,0.11764376667406434,0.11820481225140227,0.11869398247419786,0.11953867234887275,0.1205542893949228,0.12089316040520648,0.1224565487192962,0.12428970883979927,0.12675673119179087,0.12722603823517528,0.12766425078349375,0.1278281014003892,0.1283370285746668,0.12907214680385512,0.12927095215924939,0.1292699131712456,0.1290417610518484,0.1290266916584403)
sp_rate_pop_exponential_maxHPD=c(1.2436980713332384, 0.1699572178956968,0.16978317683459448,0.1699168263258951,0.1699159958933459,0.16981322361381365,0.16981322361381365,0.16994734082043716,0.16972711852543526,0.16934894771552197,0.1692983497522768,0.16930882655929413,0.16930882655929413,0.1694250253812651,0.16972677077575044,0.16972711852543526,0.1695316688120376,0.16958161572607636,0.1695316688120376,0.1695460382531998,0.16948402783200142,0.16985034093169635,0.16994734082043716,0.1700961866417923,0.1704263767022968,0.17058575531898718,0.1710285113996579,0.17137973962266848,0.17226888543111418,0.17272885303454363,0.17322962238361467,0.17389826794295077,0.17533402793631564,0.17622668410301678,0.1777442609135174,0.17975593999704825,0.1815410468091861,0.18533411308017944,0.1886925164964272,0.19052179459428867,0.19152884809195414,0.19275409929146606,0.19346607621568018,0.19475816458891332,0.19663858886221255,0.19849816466962636,0.1994522607672217,0.2003971488308411,0.2003982851607836)
#Frequency of shifts
sp_rate_pop_exponential_mids=c(-48.984375, -47.953125,-46.921875,-45.890625,-44.859375,-43.828125,-42.796875,-41.765625,-40.734375,-39.703125,-38.671875,-37.640625,-36.609375,-35.578125,-34.546875,-33.515625,-32.484375,-31.453125,-30.421875,-29.390625,-28.359375,-27.328125,-26.296875,-25.265625,-24.234375,-23.203125,-22.171875,-21.140625,-20.109375,-19.078125,-18.046875,-17.015625,-15.984375,-14.953125,-13.921875,-12.890625,-11.859375,-10.828125,-9.796875,-8.765625,-7.734375,-6.703125,-5.671875,-4.640625,-3.609375,-2.578125,-1.546875,-0.515625)
sp_rate_pop_exponential_counts=c(0.1, 0.900025,0.0014375,0.0020325,0.0019575,0.0021625,0.002075,0.00281,0.0038025,0.0022375,0.0021875,0.00233,0.00272,0.0035975,0.0034925,0.002735,0.0026275,0.0024975,0.0025775,0.0032575,0.0042875,0.00664,0.0088175,0.0069275,0.007035,0.007435,0.0074275,0.0090375,0.0092125,0.0133375,0.01397,0.0199325,0.016505,0.0218625,0.0210775,0.0211475,0.025835,0.0329475,0.03188,0.018775,0.0168325,0.0104625,0.011045,0.01177,0.015055,0.00638,0.0058275,0.0001825)
bf2 = 0.05520331253624976
bf6 = 0.3015459891543565
#Net Rate
net_rate_pop_exponential_time=c(-48.994897959183675, -47.984693877551024,-46.974489795918366,-45.964285714285715,-44.954081632653065,-43.94387755102041,-42.933673469387756,-41.923469387755105,-40.913265306122454,-39.9030612244898,-38.892857142857146,-37.88265306122449,-36.87244897959184,-35.86224489795919,-34.852040816326536,-33.84183673469388,-32.83163265306122,-31.821428571428573,-30.81122448979592,-29.801020408163268,-28.790816326530614,-27.78061224489796,-26.77040816326531,-25.760204081632654,-24.75,-23.73979591836735,-22.729591836734695,-21.71938775510204,-20.70918367346939,-19.698979591836736,-18.688775510204085,-17.67857142857143,-16.66836734693878,-15.658163265306122,-14.647959183673471,-13.63775510204082,-12.627551020408163,-11.617346938775512,-10.607142857142861,-9.59693877551021,-8.586734693877553,-7.576530612244902,-6.566326530612251,-5.556122448979593,-4.5459183673469425,-3.535714285714292,-2.525510204081634,-1.5153061224489832,-0.5051020408163254)
net_rate_pop_exponential_net_rate=c(0.8481152931071186, 0.009003756625930854,0.008956774553255524,0.008949025786747477,0.008989241202080913,0.009021507515148885,0.009102000835039946,0.009170932061588278,0.009346254478645908,0.009389512532597523,0.009447006880809721,0.009501818907991214,0.009582232525156328,0.009679294515522373,0.009802205235616715,0.009882572899794412,0.0099399118060093,0.00999658589308562,0.010048463766132335,0.010118560068510444,0.010218895984761996,0.010370959384071282,0.010666078517471709,0.010950229591983992,0.011166156430346559,0.011400534887434546,0.011622584126120753,0.011846939281733594,0.012146895049549328,0.012500638485345678,0.013012936066047462,0.013604776639611009,0.014318438919951111,0.014970707420378619,0.015785133676918996,0.01660217113238529,0.017530466761423627,0.018734034678189595,0.020421154389166987,0.02170804326506759,0.022488389425482397,0.023091051692260813,0.023469047862962392,0.023877268992724795,0.024444988914446443,0.02518748315356013,0.025416682171068777,0.025623141914423237,0.025624805400043055)
net_rate_pop_exponential_net_minHPD=c(-0.010993748875713222, -0.02003419969034638,-0.019111516609954665,-0.01875903173696221,-0.018162510702617846,-0.01735569091004266,-0.01727657932966306,-0.01727657932966306,-0.016777852395159526,-0.016330478957813754,-0.016291999518203204,-0.016116228671997662,-0.015888024335414266,-0.01564002058231001,-0.015321052605152305,-0.015226430704811972,-0.015159607792237129,-0.015147517616040718,-0.014837049407171246,-0.014768536787039196,-0.014768536787039196,-0.014350113583599347,-0.014078275101419513,-0.013609966370338414,-0.013140508756528807,-0.01274262445140345,-0.01301161952515792,-0.012795573400130184,-0.012749793030548129,-0.012700852696490572,-0.012749793030548129,-0.01230330791099124,-0.011578935459860162,-0.011337907507774636,-0.010599712091150987,-0.009613125037928383,-0.008981104721105554,-0.007977629928447333,-0.006815934892578784,-0.005014518942984075,-0.004874512943092102,-0.004083940451059143,-0.003990581034743285,-0.00432927848044154,-0.0040845044734354186,-0.0040867810180785225,-0.004773640656155792,-0.005041318832924718,-0.005057535153603465)
net_rate_pop_exponential_net_maxHPD=c(1.108603475565949, 0.03426051513576017,0.03397902830128993,0.03340886771439572,0.03323790767373293,0.03325579160537678,0.03281969919737085,0.03238737032129385,0.03208043681452495,0.032184514636400235,0.0319427115314978,0.03191412649184652,0.0319427115314978,0.03191393937009043,0.03194853339195691,0.03194853339195691,0.031946887563614634,0.03191393937009043,0.03221024493273256,0.032271498700262635,0.032271498700262635,0.03268117613806247,0.032813696566970724,0.03323868969619412,0.03374451101106821,0.034258117548342626,0.03427469170719999,0.03482808073646523,0.035374765440429196,0.036167876637918484,0.037106627047094054,0.03862923887647762,0.04055968286916903,0.04212192098570165,0.044267270872342585,0.046376392321698906,0.048539324048999216,0.0516275502656692,0.055424162532869886,0.05902552111531764,0.06024679680453571,0.062293700349414616,0.0636238357320095,0.06450955131756805,0.06647022409342093,0.068560001335228,0.06952757004123378,0.07197818697249923,0.071985487322917)
#Number Shifts
ex_rate_pop_exponential_unique=c(1.0, 2.0,3.0,4.0,5.0)
ex_rate_pop_exponential_counts=c(329810, 63555,6235,374,26)
#ex_rate_ Plot
ex_rate_pop_exponential_time=c(-48.994897959183675, -47.984693877551024,-46.974489795918366,-45.964285714285715,-44.954081632653065,-43.94387755102041,-42.933673469387756,-41.923469387755105,-40.913265306122454,-39.9030612244898,-38.892857142857146,-37.88265306122449,-36.87244897959184,-35.86224489795919,-34.852040816326536,-33.84183673469388,-32.83163265306122,-31.821428571428573,-30.81122448979592,-29.801020408163268,-28.790816326530614,-27.78061224489796,-26.77040816326531,-25.760204081632654,-24.75,-23.73979591836735,-22.729591836734695,-21.71938775510204,-20.70918367346939,-19.698979591836736,-18.688775510204085,-17.67857142857143,-16.66836734693878,-15.658163265306122,-14.647959183673471,-13.63775510204082,-12.627551020408163,-11.617346938775512,-10.607142857142861,-9.59693877551021,-8.586734693877553,-7.576530612244902,-6.566326530612251,-5.556122448979593,-4.5459183673469425,-3.535714285714292,-2.525510204081634,-1.5153061224489832,-0.5051020408163254)
ex_rate_pop_exponential_rate=c(0.13262017801660556, 0.13294980069602566,0.13299169604767042,0.13301957669095937,0.13300381859594493,0.13302919796774945,0.13300641701367616,0.1329962435766119,0.1329923251246139,0.1330022268579817,0.13300884467112925,0.1330012929992617,0.13299961371641975,0.13299948038709564,0.1330021766533782,0.13300132827818217,0.13300782629954686,0.1330138880901555,0.133017505123017,0.13302474973817338,0.1330311878296809,0.13303895798837734,0.1330473583546785,0.13305900895535902,0.13306431538737984,0.13307484931496155,0.1330931693945475,0.13311167729284612,0.13312841033484682,0.13313684296913772,0.13312912757310683,0.13313255817037617,0.13314643776186688,0.13317301077354796,0.1332067010012222,0.1332299605620627,0.13325470144370657,0.13328287864610647,0.1333246024259173,0.13339000105925772,0.13346841767889306,0.13358621841134774,0.13371346757767888,0.1338355359137384,0.13401599110350737,0.13416538672698222,0.13424106944278572,0.13435594820762084,0.13435704806007212)
ex_rate_pop_exponential_minHPD=c(0.11146264852517329, 0.11274457965468654,0.11345960745047577,0.11332195447007196,0.11380572880862028,0.11368782293171124,0.11397361327416507,0.11423135793651616,0.11448027079769751,0.11467706640419903,0.11466957242040632,0.11462503261791314,0.11475464263790837,0.11466026852655059,0.11475191610362745,0.1147070175531622,0.11470114343673762,0.1147513568389307,0.11479240162360542,0.11485073110256706,0.11487163229237611,0.11485044450309484,0.11487163229237611,0.11487173057158813,0.11510483056354073,0.11486966317436266,0.11487163229237611,0.11509977615695746,0.11509151059544305,0.11509144833913292,0.1149498488666781,0.11510587538089136,0.11508989873753989,0.11487130559214292,0.1149502544528071,0.11492635436512538,0.11487146282126423,0.11487130559214292,0.11492635436512538,0.11481767667958144,0.11480562917104875,0.11502633717751207,0.11487146282126423,0.11462021985124882,0.11460812326348233,0.11436050278778852,0.11452627723081699,0.11393943210571642,0.11391853411438052)
ex_rate_pop_exponential_maxHPD=c(0.15764863547696395, 0.1553554047984202,0.15523168912738824,0.1544396072719023,0.15439341982960128,0.15377981440514613,0.1537829857955684,0.15379279933355358,0.15377676583583708,0.15377676583583708,0.15361079045321124,0.15341444462948967,0.15339983320707176,0.15320000954111315,0.15319249818165276,0.15309412354967938,0.15303960490128843,0.15304020114859934,0.15303856084818523,0.15303960490128843,0.15300882897793197,0.15295011215078008,0.15294686200487426,0.15293080622253497,0.1531689227025785,0.15294792020257608,0.15294758166697953,0.1531689227025785,0.15316761259920791,0.15323075940338288,0.15320048355856547,0.15341907157701482,0.1535009783854792,0.15341907157701482,0.1536065526336192,0.15361109774149434,0.15358923645862738,0.15364128451559186,0.15377676583583708,0.15377723085643974,0.15399094277287131,0.1544593635513917,0.154558632862162,0.154677794040954,0.15526318143106224,0.15577564221004653,0.15656647406243795,0.15696249284258446,0.15695199612738678)
#Frequency of shifts
ex_rate_pop_exponential_mids=c(-48.984375, -47.953125,-46.921875,-45.890625,-44.859375,-43.828125,-42.796875,-41.765625,-40.734375,-39.703125,-38.671875,-37.640625,-36.609375,-35.578125,-34.546875,-33.515625,-32.484375,-31.453125,-30.421875,-29.390625,-28.359375,-27.328125,-26.296875,-25.265625,-24.234375,-23.203125,-22.171875,-21.140625,-20.109375,-19.078125,-18.046875,-17.015625,-15.984375,-14.953125,-13.921875,-12.890625,-11.859375,-10.828125,-9.796875,-8.765625,-7.734375,-6.703125,-5.671875,-4.640625,-3.609375,-2.578125,-1.546875,-0.515625)
ex_rate_pop_exponential_counts=c(0.0004475, 0.0159025,0.0063825,0.00588,0.0053475,0.005125,0.00383,0.0037475,0.0036725,0.0036775,0.003255,0.003265,0.003425,0.0033075,0.003225,0.002835,0.002615,0.0026075,0.00264,0.002855,0.0027075,0.0025625,0.002525,0.002615,0.002825,0.00287,0.0028175,0.0030375,0.003015,0.0036875,0.0029625,0.0031675,0.00375,0.003865,0.00304,0.00297,0.0031775,0.0034875,0.0040925,0.0042825,0.005505,0.006305,0.0056625,0.0069525,0.0071525,0.00642,0.0073375,0.0002925)
#Net Diversity
diversity_pop_exponential_time=c(-48.994897959183675, -47.984693877551024,-46.974489795918366,-45.964285714285715,-44.954081632653065,-43.94387755102041,-42.933673469387756,-41.923469387755105,-40.913265306122454,-39.9030612244898,-38.892857142857146,-37.88265306122449,-36.87244897959184,-35.86224489795919,-34.852040816326536,-33.84183673469388,-32.83163265306122,-31.821428571428573,-30.81122448979592,-29.801020408163268,-28.790816326530614,-27.78061224489796,-26.77040816326531,-25.760204081632654,-24.75,-23.73979591836735,-22.729591836734695,-21.71938775510204,-20.70918367346939,-19.698979591836736,-18.688775510204085,-17.67857142857143,-16.66836734693878,-15.658163265306122,-14.647959183673471,-13.63775510204082,-12.627551020408163,-11.617346938775512,-10.607142857142861,-9.59693877551021,-8.586734693877553,-7.576530612244902,-6.566326530612251,-5.556122448979593,-4.5459183673469425,-3.535714285714292,-2.525510204081634,-1.5153061224489832,-0.5051020408163254)
diversity_pop_exponential_net_diversity=c(77.77, 77.96,78.17,78.13,77.35,77.7,77.69,77.02,78.2,78.87,78.79,78.38,78.93,79.28,79.99,79.59,79.73,80.16,80.27,80.39,80.28,80.81,81.41,82.54,82.73,83.69,84.54,85.38,85.67,86.29,88.32,89.26,90.79,92.07,93.82,96.05,97.93,100.0,102.41,104.5,109.36,112.82,117.29,121.38,126.02,130.68,136.13,141.81,148.07)

#n <- dev.off()