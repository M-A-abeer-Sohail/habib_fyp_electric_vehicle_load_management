function newVi = getNewVi(Vpu1hr, ViSetting, ViLen)
    newVi = zeros(1,ViLen);
    for i = 1:ViLen
        newVi(i) = Vpu1hr(ViSetting(i,1),ViSetting(i,2));
    end
end