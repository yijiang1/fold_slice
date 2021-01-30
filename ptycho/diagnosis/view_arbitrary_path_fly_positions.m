%view_arbitrary_path_fly_positions.m
%show scan positions in arbitrary-path fly-scan reconstruction

%% load a matlab recon file
figure
hold on
for i=1:size(outputs.probe_positions,3)
    scatter(outputs.probe_positions(:,1,i)*p.dx_spec(1),outputs.probe_positions(:,2,i)*p.dx_spec(1),'.','DisplayName',num2str(i)); axis image
end
legend