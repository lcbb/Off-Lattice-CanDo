% Check whether a finite element simulation converged.
function status = checkConvergence(msgPATH)

status = 0;

fid = fopen(msgPATH);

tline = fgetl(fid);
while(ischar(tline))
    if(strfind(tline,'ENDCODE=1'))
        status = -1;
        fclose(fid);
        return;
    end
    tline = fgetl(fid);
end

fclose(fid);

end