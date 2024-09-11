function [] = make_table(ResFolder,GE_cond,Outputs)

FID = fopen(fullfile(ResFolder,'Table1.tex'),'w');
fprintf(FID,' \\begin{tabular}{lcc} \\hline \n');
fprintf(FID,'  & US Data & Model \\\\ \n');
fprintf(FID,' \\hline \n');

fprintf(FID,'%s  & %8.3f & %8.3f \\\\ \n','Top 10 Employment',0.67,Outputs.top10_empl);
fprintf(FID,'%s  & %8.3f & %8.3f \\\\ \n','Top 5 Earnings',0.30,Outputs.top5_earnings);
fprintf(FID,'%s  & %8.3f & %8.3f \\\\ \n','Establisments exit rate',0.10,Outputs.exit_E_to_W);
fprintf(FID,'%s  & %8.3f & %8.3f \\\\ \n','Real interest rate',0.045,Outputs.r);

fprintf(FID, '\\hline \n \\end{tabular} \n');
fclose(FID);

%------ On screen
disp('---------------------------')
disp('CALIBRATION TABLE')

width = length('Establisments exit rate');

fprintf('                          US Data  Model \n');
fprintf('%-*s   %-8.4f %-8.4f \n',width,'Top 10 Employment',0.67,Outputs.top10_empl);
fprintf('%-*s   %-8.4f %-8.4f \n',width,'Top 5 Earnings',0.30,Outputs.top5_earnings);
fprintf('%-*s   %-8.4f %-8.4f \n',width,'Establisments exit rate',0.10,Outputs.exit_E_to_W);
fprintf('%-*s   %-8.4f %-8.4f \n',width,'Real interest rate',0.045,Outputs.r);

fprintf('%-*s   %-8.4f %-8.4f \n',width,'Real interest rate',0.045,Outputs.r);
fprintf('%-*s   %-8.4f %-8.4f \n',width,'GE condition 1',0.0,GE_cond(1));
fprintf('%-*s   %-8.4f %-8.4f \n',width,'GE condition 2',0.0,GE_cond(2));
disp('---------------------------')

end %end function
