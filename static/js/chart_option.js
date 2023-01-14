$(document).ready(function()
{
    $('#chart_choice').on('change','select',function(){
        var select_value = document.getElementById('chart_type').value;
        var y_axis = document.getElementById('y-div');
        if((select_value==0 || select_value==1) && !y_axis)
        {
            var axis_choice = document.getElementById('axis-choice');
            var y_s
            return;
        }
        if(select_value==2 && y_axis)
        {
            y_axis.remove();
            return;
        }

    });
});