%Edit From Harry
function PreData = ARIMA(SourceData,step)
  %warning off all
  TempData=SourceData;
  TempData=detrend(TempData);%ȥ������
  TrendData=SourceData-TempData;%���ƺ���
  %---------------------------------------------------��֣�ƽ�Ȼ�ʱ������---------
  [H,PValue,TestStat,CriticalValue] = adftest(TempData);
  %difftime=0;
  SaveDiffData=[];
  while ~H
      SaveDiffData=[SaveDiffData,TempData(1,1)];
      TempData=diff(TempData);%��֣�ƽ�Ȼ�ʱ������
      %difftime=difftime+1;%��ִ���
      [H,PValue,TestStat,CriticalValue] = adftest(TempData);%adf���飬�ж�ʱ�������Ƿ�ƽ�Ȼ�
  end
  %---------------------------------------------------ģ�Ͷ��׻�ʶ��--------------
  u = iddata(TempData');
  test = [];
  for p = 0:5    %�Իع��ӦPACF,�����ͺ󳤶�����p��q��һ��ȡΪT/10��ln(T)��T^(1/2),����ȡT/10=12
    for q = 0:5    %�ƶ�ƽ����ӦACF
        if(p + q ~= 0)  %added by Harry on 20110214 begin
            m = armax(u,[p q]);        
            AIC = aic(m);  %armax(p,q),����AIC
            test = [test;p q AIC];
        end
    end
  end
  for k = 1:size(test,1)
    if test(k,3) == min(test(:,3)) %ѡ��AICֵ��С��ģ��
        p_test = test(k,1);
        q_test = test(k,2);
        break;
    end
  end
  %---------------------------------------------------1��Ԥ��-----------------
  for j=1:step
      TempData=[TempData,0];
  end
  n=iddata(TempData');
  m = armax(u,[p_test q_test]);        %armax(p,q),[p_test q_test]��ӦAICֵ��С
  P1=predict(m,n,step);
  PreR=P1.OutputData;
  PreR=PreR';
  %---------------------------------------------------��ԭ���-----------------
  if size(SaveDiffData,2)~=0
      for index=size(SaveDiffData,2):-1:1
          PreR=cumsum([SaveDiffData(index),PreR]);
      end
  end 
  %---------------------------------------------------Ԥ�����Ʋ����ؽ��----------------
  mp1=polyfit([1:size(TrendData,2)],TrendData,1);
  xt=[];
  for j=1:step
      xt=[xt,size(TrendData,2)+j];
  end
  TrendResult=polyval(mp1,xt);
  PreData=TrendResult+PreR(size(SourceData,2)+1:size(PreR,2));
  %tempx=[TrendData,TrendResult]+PreR;
  %plot(tempx,'r'),hold on,plot(SourceData);
  

  
  
  
     