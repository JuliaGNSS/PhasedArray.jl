using Test, Statistics, LinearAlgebra, PhasedArray, StaticArrays, CoordinateTransformations, Interpolations, Random

const EXAMPLE_LUT = zeros(Complex{Float64}, 2, 10, 20)
EXAMPLE_LUT[1,:,:] = [-5.7078+9.35399im -2.5379+10.66im 0.880432+10.9225im 4.21258+10.1158im 7.13237+8.31898im 9.35399+5.7078im 10.66+2.5379im 10.9225-0.880432im 10.1158-4.21258im 8.31898-7.13237im 5.7078-9.35399im 2.5379-10.66im -0.880432-10.9225im -4.21258-10.1158im -7.13237-8.31898im -9.35399-5.7078im -10.66-2.5379im -10.9225+0.880432im -10.1158+4.21258im -8.31898+7.13237im; -2.92983+11.3667im 1.89319+11.5399im 5.95279+9.87774im 8.67235+7.17072im 10.0222+4.2519im 10.3557+1.69283im 10.1422-0.309334im 9.74-1.86734im 9.2859-3.24392im 8.70125-4.68534im 7.77128-6.30552im 6.24428-8.03485im 3.92146-9.61231im 0.748499-10.6018im -3.07119-10.4589im -6.97975-8.70517im -10.101-5.20192im -11.5001-0.382397im -10.5936+4.75934im -7.47524+8.99718im; -0.216611+11.5061im 5.35705+10.2762im 9.10426+6.86225im 10.5977+2.85826im 10.3932-0.408505im 9.46519-2.30716im 8.64987-2.92284im 8.33038-2.79908im 8.44385-2.58407im 8.67441-2.7732im 8.65067-3.60217im 8.05841-5.08176im 6.62685-7.05999im 4.08195-9.14687im 0.275511-10.5469im -4.41347-10.1469im -8.77984-7.11999im -11.1258-1.69957im -10.2009+4.57227im -6.07889+9.51627im;1.70403+9.96745im 6.80512+7.7488im 9.4351+3.47579im 9.59045-0.925037im 8.37188-4.06076im 7.08296-5.34748im 6.55379-5.01314im 6.84712-3.83431im 7.54192-2.63336im 8.1795-1.97996im 8.42301-2.10676im 8.09846-3.03599im 7.0956-4.77848im 5.1519-7.18815im 1.87367-9.50075im -2.66495-10.2227im -7.25878-7.95956im -9.83949-2.67411im -8.80547+3.7818im -4.31549+8.63307im; 2.58851+7.52972im 6.36233+5.13796im 7.60694+1.07401im 6.71037-2.98171im 5.09368-5.87641im 4.10786-6.99832im 4.31087-6.42438im 5.32356-4.88562im 6.4295-3.26498im 7.2145-2.17533im 7.46964-1.83023im 7.10368-2.23904im 6.1511-3.53087im 4.49189-5.75432im 1.76965-8.25608im -2.03981-9.48682im -5.93235-7.862im -7.97165-3.20938im -6.72227+2.58651im -2.53476+6.74812im; 2.65306+4.98345im 4.89296+3.14456im 4.94449+0.00609865im 3.40884-3.31952im 1.78895-6.01659im 1.32315-7.3719im 2.21557-7.11829im 3.74402-5.72856im 5.0443-4.08801im 5.81967-2.86509im 6.02191-2.27207im 5.56169-2.31398im 4.48553-3.1823im 2.8431-5.00265im 0.546714-7.21154im -2.28902-8.41962im -4.9136-7.24029im -5.93083-3.53102im -4.46613+1.10049im -1.04682+4.38743im; 2.28216+2.76921im 3.19345+1.88159im 2.43028-0.0756372im 0.656211-2.57115im -0.778546-5.07559im -0.835239-6.77379im 0.444139-7.05616im 2.17514-6.0952im 3.48497-4.67232im 4.13565-3.48093im 4.19436-2.78641im 3.61634-2.62313im 2.41185-3.16898im 0.750127-4.5462im -1.11641-6.28675im -2.89148-7.2718im -4.08302-6.49968im -4.01169-3.90781im -2.44014-0.572563im 0.0448655+1.97723im; 1.83884+1.02319im 1.77609+1.1203im 0.590593+0.257649im -1.09378-1.43161im -2.24822-3.58215im -2.11921-5.39966im -0.838565-6.13043im 0.763675-5.69152im 1.92109-4.62517im 2.38188-3.54888im 2.20869-2.81842im 1.46806-2.56576im 0.233282-2.91209im -1.24983-3.91966im -2.53723-5.24068im -3.26216-6.09047im -3.22145-5.81225im -2.34342-4.32608im -0.824032-2.15118im 0.801839-0.134057im; 1.55077-0.169921im 0.990389+0.583706im -0.223209+0.520155im -1.57787-0.399746im -2.41894-1.90247im -2.32723-3.38619im -1.42979-4.2352im -0.288417-4.2298im 0.541021-3.60663im 0.791734-2.7882im 0.456765-2.14186im -0.337265-1.89663im -1.39584-2.16497im -2.4216-2.93952im -3.03034-3.97791im -2.95759-4.80743im -2.2056-4.99035im -0.998067-4.36148im 0.304165-3.07186im 1.27485-1.51534im; 1.45956-0.605368im 0.958095+0.206997im 0.123812+0.541103im -0.738487+0.333007im -1.32277-0.27992im -1.46311-1.00807im -1.21287-1.54739im -0.799262-1.72811im -0.487362-1.57483im -0.455091-1.25732im -0.741344-0.993681im -1.25664-0.965799im -1.82185-1.26631im -2.21663-1.86893im -2.24503-2.61925im -1.81619-3.27032im -0.996652-3.57269im 0.00473669-3.37546im 0.906467-2.68569im 1.44455-1.67001im]

EXAMPLE_LUT[2,:,:] = [-6.05477+9.15928im -2.92805+10.582im 0.485278+10.9689im 3.85111+10.2821im 6.83996+8.5888im 9.15928+6.05477im 10.582+2.92805im 10.9689-0.485278im 10.2821-3.85111im 8.5888-6.83996im 6.05477-9.15928im 2.92805-10.582im -0.485278-10.9689im -3.85111-10.2821im -6.83996-8.5888im -9.15928-6.05477im -10.582-2.92805im -10.9689+0.485278im -10.2821+3.85111im -8.5888+6.83996im; -1.85229+10.2583im 0.125939+10.0823im 1.65825+9.72805im 3.02163+9.33267im 4.46866+8.81454im 6.11648+7.95984im 7.9025+6.51552im 9.57352+4.27158im 10.6959+1.1526im 10.7169-2.66578im 9.13052-6.65172im 5.75145-9.93518im 0.963961-11.551im -4.25862-10.8537im -8.66404-7.87352im -11.2267-3.36423im -11.5572+1.50995im -9.98569+5.66294im -7.31158+8.47201im -4.40008+9.88433im; 2.32966+9.42276im 2.90401+8.58528im 2.69771+8.27122im 2.41763+8.40792im 2.57174+8.66968im 3.39564+8.71228im 4.91768+8.22935im 6.98437+6.92317im 9.19164+4.49905im 10.7563+0.780573im 10.5653-3.91931im 7.72874-8.44848im 2.38227-11.0803im -3.99127-10.4541im -9.1787-6.52537im -11.4405-0.695218im -10.3918+4.98253im -7.01941+8.89024im -2.95495+10.5135im 0.395107+10.3614im; 5.49062+7.08307im 5.11083+6.5063im 3.78307+6.77099im 2.44312+7.46082im 1.70405+8.09974im 1.81853+8.40605im 2.83353+8.21762im 4.71836+7.37907im 7.26585+5.60874im 9.73579+2.493im 10.6645-1.99877im 8.5916-6.75894im 3.37442-9.68142im -3.19291-9.01066im -8.30866-4.77566im -9.94985+1.1705im -7.91714+6.37193im -3.64795+9.18597im 0.860986+9.5092im 4.12675+8.37315im; 7.22412+4.17167im 6.57231+4.34252im 4.83524+5.29677im 3.00893+6.36361im 1.79323+7.10101im 1.45678+7.36378im 2.01696+7.11406im 3.51693+6.36654im 5.92561+4.94878im 8.59626+2.45252im 10.0008-1.26968im 8.51171-5.33284im 3.86572-7.75216im -2.0734-6.89915im -6.4878-3.00572im -7.5591+2.00388im -5.31414+5.86753im -1.19521+7.29844im 3.00733+6.59293im 6.04548+5.10292im; 7.68445+1.47502im 7.27835+2.39481im 5.64459+3.86568im 3.74469+5.08031im 2.3759+5.73246im 1.82336+5.83397im 2.08756+5.41418im 3.24326+4.55491im 5.30984+3.213im 7.69939+1.18895im 9.03165-1.53379im 7.90794-4.31362im 4.1157-5.69674im -0.712463-4.62099im -4.23839-1.5051im -5.06215+2.04657im -3.29416+4.34303im -0.0300219+4.5657im 3.49157+3.23854im 6.33615+1.80823im; 7.19303-0.570016im 7.24467+0.814603im 5.97234+2.51897im 4.23478+3.68909im 2.87922+4.12406im 2.24873+3.96792im 2.36568+3.31865im 3.29102+2.27606im 4.98074+0.930567im 6.90681-0.636702im 7.97693-2.27709im 7.17668-3.57541im 4.41031-3.79126im 0.830211-2.55616im -1.93561-0.37611im -2.8617+1.67763im -1.95045+2.59107im 0.197472+1.9669im 2.93443+0.426366im 5.56805-0.733088im; 5.91623-1.7095im 6.34813-0.259612im 5.52823+1.34498im 4.09749+2.32072im 2.83516+2.47207im 2.18118+1.98769im 2.24191+1.0704im 3.02492-0.106705im 4.39591-1.33057im 5.91086-2.33402im 6.81313-2.90305im 6.44273-2.88759im 4.73381-2.17793im 2.30921-0.910476im 0.129227+0.432054im -1.05343+1.26404im -1.03421+1.1553im 0.0635849+0.0774827im 2.00228-1.34476im 4.23862-2.14238im; 3.9379-1.77683im 4.44186-0.675618im 4.0173+0.476719im 3.01398+1.10376im 1.9985+0.992059im 1.42103+0.27258im 1.47752-0.786439im 2.1598-1.90524im 3.29831-2.79855im 4.54519-3.19518im 5.41629-2.95136im 5.50112-2.13123im 4.6906-0.968832im 3.22474+0.19088im 1.58417+0.96663im 0.286035+1.08017im -0.305234+0.466316im -0.0127989-0.651365im 1.10377-1.74274im 2.63825-2.21163im; 1.45027-0.844732im 1.68234-0.438207im 1.50331-0.0300361im 1.03215+0.101834im 0.530564-0.177733im 0.268351-0.817607im 0.413718-1.62702im 0.988208-2.35709im 1.87145-2.77413im 2.83284-2.72152im 3.59225-2.17278im 3.91054-1.25277im 3.67653-0.204476im 2.94852+0.686521im 1.93736+1.18077im 0.940967+1.16338im 0.249812+0.689672im 0.0478472-0.0214697im 0.334226-0.661207im 0.907742-0.964745im]

Random.seed!(1234)

include("manifold.jl")
include("plots.jl")
include("filter.jl")
