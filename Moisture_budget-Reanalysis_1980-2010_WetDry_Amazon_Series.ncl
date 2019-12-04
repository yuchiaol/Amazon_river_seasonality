;=======================================================================

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;=======================================================================
; Output mrso data and resample (soil moisture content)

begin

 dnum     = 6 
 varnum   = 11 
 PI       = 3.14159265359
 nx       = 360
 ny       = 180
 nyear    = 32        

 var        = (/"Precip","Precip","SLHF","vdq","wdq","res","wbdqp","wpdqb","wpdqp","W500","Qint","wbdqb"/)
 prname     = (/"GPCP","GPCC","PRECL","CMAP","CRU","MSWEP"/)
 oendyear   = (/2018,2016,2018,2014,2016,2014/)

 ;======================================================================
 ; Read File year star and end (ifsy ifey ofsy ofey)

 dataname = new(dnum,"string")
 end_year = new(dnum,"integer")

 pwd1 = "/home/L.JeffBlue/work2/62.Amazon_Streamflow/5th-Moisture_Budget/data/"
 inputfilename = pwd1+"/Data_file_year_res.dat"

 data  = asciiread(inputfilename, -1, "string")
 cdata = stringtochar(data)
 dataname = str_right_strip(chartostring(cdata(0:dnum-1,0:19)))

 dd   = chartostring(cdata(0:dnum-1,25:28))
 end_year = stringtointeger(dd)

 delete([/data,cdata,dd/])

 print(""+dataname)

 ;======================================================================

 rcva   = new((/3,varnum/),"float",1e+20)    ;WetDryRange,Model,GlobalOceanLand,Var
 rcvi   = new((/3,varnum/),"float",1e+20)
 proba  = new((/3,varnum/),"float",1e+20)
 probj  = new((/3,varnum/),"float",1e+20)

 ts_dd    = new((/3,varnum,dnum+1,nyear-1/),"float",1e+20)
 trend_dd = new((/3,varnum,nyear-1/),"float",1e+20)
 dd       = new((/3,nyear-1/),"float",1e+20)
 unit     = new(varnum,"string")

 dnum = 3

 do v=0,0
   do d=0,dnum-1
     print("========"+prname(d)+" "+var(v)+"========")
     print("======Read Data======")

     pwd1 = "/home/L.JeffBlue/work2/62.Amazon_Streamflow/3rd-OBSPR/data/"

     ifn1 = pwd1+"1-"+prname(d)+"_1979-"+oendyear(d)+"_Precip_WetDry_Season_RegAVG_Data_3rmon.bin"

     scat = 28.0

     dyear = oendyear(d)-1979+1

     tsdd = fbindirread(ifn1,0,(/3,dyear/),"float")

     do s=0,2
       ts_dd(s,v,d,:) = tsdd(s,1:nyear-1)*scat
     end do

     delete(tsdd)
   end do
 end do

 dnum = 6

 do v=1,varnum-1
   ;=======================================================================
   ; read data

   do d=0,dnum-1
     print("========"+dataname(d)+" "+var(v)+"========")
     print("======Read Data======")

     pwd1 = "/home/L.JeffBlue/work2/62.Amazon_Streamflow/5th-Moisture_Budget/data/"

     if (v.eq.1) then
       ifn1 = pwd1+"9-"+dataname(d)+"_1979-2010_"+var(v)+"_3rm_WetDryRan_Rainforest_Avg.dat"
       scat = 28.0
     else if (v.eq.9) then
       ifn1 = pwd1+"9-"+dataname(d)+"_1979-2010_"+var(v)+"_3rm_WetDryRan_Rainforest_Avg.dat"
       scat = 100.0
     else
       ifn1 = pwd1+"9-"+dataname(d)+"_1979-2010_"+var(v)+"_3rm_WetDryRan_Rainforest_Avg.dat"
       scat = 1.0
     end if
     end if

     data = asciiread(ifn1,(/nyear,6/),"float")
     data = data*scat

     do s=0,2
       ts_dd(s,v,d,:) = data(1:nyear-1,s)
     end do
   end do
 end do

 ts_dd(:,:,dnum,:) = dim_avg_n(ts_dd,2)

   print("======Linear Regression Coeficient======")

 do v=0,varnum-1
   do s=0,2
     do d=0,dnum
       x = fspan(1,nyear-1,nyear-1)
       y = ts_dd(s,v,d,:)

       rc = regline_stats(x/10,y)

       ;trend_dd(s,v,0,:) = rc@Yest
       ;trend_dd(s,v,1,:) = rc@YMR025
       ;trend_dd(s,v,2,:) = rc@YMR975
       ;trend_dd(s,v,1,:) = rc@YPI025
       ;trend_dd(s,v,2,:) = rc@YPI975
       ;trend_dd(s,v,1,:) = rc@b95(0)*x/10 + rc@y95(0)
       ;trend_dd(s,v,2,:) = rc@b95(1)*x/10 + rc@y95(1)

       proba(s,v) = rc@pval(1)
       rcva(s,v)  = rc@b(1)
       rcvi(s,v)  = (rc@b95(1)-rc@b95(0))/2.0

       trend_dd(s,v,:) = rc@b(1)*x/10+rc@b(0)

       ;df         = rc@nptxy-2
       ;proba(s,v) = betainc(df/(df+rc@tval(1)^2), df/2.0, 0.5)

       ;============================P value adjust======================

       N          = rc@nptxy
       sigr       = 0.05
       df         = equiv_sample_size (y, sigr,0) 
       probj(s,v) = betainc(df/(df+rc@tval(1)^2), df/2.0, 0.5)

       aa = 0 
       if (.not.ismissing(proba(s,v)).AND.proba(s,v).gt.0.05) then
         aa = 1
       end if
       bb = 0
       if (.not.ismissing(probj(s,v)).AND.probj(s,v).gt.0.05) then
         bb = 1
       end if

       print("=="+d+"== "+var(v)+" "+s+" Rc="+rcva(s,v)+" Inter="+rc@yintercept+"  "+proba(s,v))
       if (aa.ne.bb) then
         print("=="+d+"== "+var(v)+" "+s+" "+proba(s,v)+" "+probj(s,v)+" "+df+" "+N+"==")
       end if
     end do
     print("==================")
   end do
 end do

 varname = (/"Precipitation","Precipitation","Evapotranspiration","Horizontal Advection",\
             "Vertical Advection","Residual","Thermodynamical Component","Dynamical Component",\
             "Non-linear Component","Vertical Velocity at 500hPa","Total Column Water Vapor",\
             "-<~F33~w~F21~~H-24V12F18~-~FH0V-12~~F34~6~F10~~B~p~N~q~H-24V12F18~-~FH0V-12~~F21~>"/)

 subfig   = (/"a","a","b","c","d","e","f","g","h","i","j","k","l"/)
 unit     = "W/m~S~2~N~"
 unit(9)  = "100*Pa/s"
 unit(10) = "Kg/Kg" 
 season   = (/"Wet","Dry","Annual_Range"/)

 xaix = ispan(1980,2010,1)

 ;====================================================================
 ; Create Plot

 do s=0,2
   plot  = new(varnum,"graphic")
   pindex = s+1

   figname = "Fig1."+pindex+"-Reanalysis_1980-2010_"+season(s)+"_Amazon_Series"
   figtype = "oldps"

   wks  = gsn_open_wks (figtype, figname)
   setvalues wks
     "wkForegroundColor" : (/0.,0.,0./)
     "wkBackgroundColor" : (/1.,1.,1./)
   end setvalues

 ;======================================================================
 ; Resource
 
   res                              = True                         ; plot mods desired

   res@gsnDraw                      = False                       ; don't draw yet
   res@gsnFrame                     = False                       ; don't advance frame yet

   res@vpWidthF                     = 0.8
   res@vpHeightF                    = 0.2

   res@tmXTOn                       = False  ; no XR tick marks
   res@tmYROn                       = True  ; no YR tick marks

   res@trXMaxF                      = 2010                          ; max value on y-axis
   res@trXMinF                      = 1980
   res@trYMaxF                      =  0.0                        ; max value on y-axis
   res@trYMinF                      =  0.0

   ;res@tmXBMode                     = "Explicit"         ; Define own tick mark labels
   ;res@tmXBValues                   = xtik
   ;res@tmXBLabels                   = xlab
   res@tmXBMajorOutwardLengthF      = -0.0
   res@tmXBMinorOutwardLengthF      = -0.0
   res@tmYLMajorOutwardLengthF      = -0.0
   res@tmYLMinorOutwardLengthF      = -0.0
   res@tmXBLabelFont                = 21                                 ; tick mark setting
   res@tmYLLabelFont                = 21
   res@tmYRLabelFont                = 21
   res@tmXBLabelFontHeightF         = 0.02
   res@tmYLLabelFontHeightF         = 0.02
   res@tmYRLabelFontHeightF         = 0.02

   res@tiXAxisString                = " "    ; xaxis
   res@tiYAxisString                = " "                  ; yaxis
   res@tiMainFont                   = 21
   res@tiXAxisFont                  = 21
   res@tiYAxisFont                  = 21
   res@tiMainFontHeightF            = 0.015
   res@tiXAxisFontHeightF           = 0.015
   res@tiYAxisFontHeightF           = 0.015

   res@gsnLeftString                = " "
   res@gsnLeftStringFontHeightF     = 0.02
   res@gsnCenterString              = " " ; add the gsn titles
   res@gsnCenterStringFontHeightF   = 0.025
   res@gsnRightString               = " "
   res@gsnRightStringFontHeightF    = 0.015
   res@gsnStringFont                = 21

 ;======================================================================

   do v=0,varnum-1
     res@xyLineThicknesses            = (/2.2,1,1/)                   ; Define line thicknesses
     res@xyDashPatterns               = (/0,0,0/)                   ; Dash patterns
     res@xyLineColors                 = (/"black","lightskyblue1","lightskyblue1"/)

     res@gsnXYAboveFillColors         = "lightskyblue1"
     res@gsnXYBelowFillColors         = "lightskyblue1"

     res@trYMaxF                      = max(ts_dd(s,v,dnum,:)+dim_stddev_n(ts_dd(s,v,0:dnum-1,:),0))
     res@trYMinF                      = min(ts_dd(s,v,dnum,:)-dim_stddev_n(ts_dd(s,v,0:dnum-1,:),0))

     res@gsnLeftString                = "("+subfig(v)+") "+varname(v)

     perc_trend = rcva(s,v)/rcva(s,0)*100

     if (v.lt.2.OR.v.gt.8) then
       if (proba(s,v).le.0.05) then
         res@gsnRightString      = "*Trend= "+sprintf("%4.2f",rcva(s,v))+"~F34~1~F21~"+sprintf("%4.2f",rcvi(s,v))+\
                                   "("+unit(v)+"/10y)"
       else
         res@gsnRightString      = "Trend= "+sprintf("%4.2f",rcva(s,v))+"~F34~1~F21~"+sprintf("%4.2f",rcvi(s,v))+\
                                   "("+unit(v)+"/10y)"
       end if
     end if

     if (v.ge.2.AND.v.le.8) then
       if (proba(s,v).le.0.05) then
         res@gsnRightString      = "*Trend= "+sprintf("%4.2f",rcva(s,v))+"~F34~1~F21~"+sprintf("%4.2f",rcvi(s,v))+\
                                   "("+unit(v)+"/10y)("+sprintf("%5.1f",perc_trend)+"%)"
       else
         res@gsnRightString      = "Trend= "+sprintf("%4.2f",rcva(s,v))+"~F34~1~F21~"+sprintf("%4.2f",rcvi(s,v))+\
                                   "("+unit(v)+"/10y)("+sprintf("%5.1f",perc_trend)+"%)"
       end if
     end if

     dd(0,:) = ts_dd(s,v,dnum,:)
     dd(1,:) = ts_dd(s,v,dnum,:)+dim_stddev_n(ts_dd(s,v,0:dnum-1,:),0)
     dd(2,:) = ts_dd(s,v,dnum,:)-dim_stddev_n(ts_dd(s,v,0:dnum-1,:),0)

     plot(v) = gsn_csm_xy(wks,xaix,dd,res)

     res@xyLineThicknesses            = (/1.5,1,1/)                   ; Define line thicknesses
     res@xyDashPatterns               = (/2,0,0/)                   ; Dash patterns
     res@xyLineColors                 = (/"red","black","grey70"/)

     delete([/res@gsnXYAboveFillColors,res@gsnXYBelowFillColors/])

     dd(0,:) = trend_dd(s,v,:)

     plot1 = gsn_csm_xy(wks,xaix,dd(0,:),res)

     overlay(plot(v),plot1)
   end do
 
 ;======================================================================
 ; Create Panel Plot

   resP                      = True                                      ; modify the panel plot
   resP@gsnFrame             = False                                     ; don't advance panel plot

   resP@txString             = "~F22~1980-2010 Reanalysis Moisture Budget Series ("+season(s)+")"
   resP@txFontHeightF        = 0.02

   gsn_panel(wks,plot(1:10),(/5,2/),resP)                               ; now draw as one plot

   frame(wks)

 ;========================Transform Image File========================

   if (figtype.eq."oldps") then
       figtype = "jpg"
       system("convert -density 800 -trim "+figname+".ps  "+figname+"."+figtype)
       system("/bin/rm "+figname+".ps")
   end if
   print("Output Figure Name = "+figname+"."+figtype)

 ;====================================================================
 end do
end

