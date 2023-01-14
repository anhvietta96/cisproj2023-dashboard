/*ApexCharts Lib, functional but deprecated */
/*
var dynamicchart_bool = document.getElementById("dynamic-chart");
if(dynamicchart_bool){
  const options = JSON.parse(document.currentScript.nextElementSibling.textContent);
  console.log(options);
  var chart = new ApexCharts(document.querySelector('#dynamic-chart'),options);
  chart.render();
}
*/

function checkValidInt(input)
{
 return typeof(input)=="string" && !isNaN(input) && !isNaN(parseInt(input));
}

function checkValidFloat(input)
{
  return typeof(input)=="string" && !isNaN(input) && !isNaN(parseFloat(input));
}

function numerical_filter(data,pos,lower_bound,upper_bound)
{
  var filtered_data = [];
  for(const i of Object.keys(data))
  {
    var sub_arr = [];
    for(let j = 0;j<data[i].length;j++)
    {
      if(data[i][j][pos]>lower_bound && data[i][j][pos]<upper_bound)
      {
        sub_arr.push(data[i][j]);
      }
    }
    filtered_data.push(sub_arr);
  }
  return filtered_data;
}

function lexicographic_filter(data,pos,substring)
{
  var filtered_data = [];
  for(const i of Object.keys(data))
  {
    var sub_arr = [];
    for(let j = 0;j<data[i].length;j++)
    {
      if(data[i][j][pos].includes(substring))
      {
        sub_arr.push(data[i][j]);
      }
    }
    filtered_data.push(sub_arr);
  }
  return filtered_data;
}

function get_selected_inchikey(option,chart_type)
{
  if([1,2].includes(chart_type))
  {
    var position = 2;
  }
  else if([3].includes(chart_type))
  {
    var position = 1;
  }
  var inchikey_array = [];
  for(const group of Object.keys(option['series']))
  {
    for(let i = 0; i < option['series'][group]['data'].length; i++)
    {
      inchikey_array.push(option['series'][group]['data'][i][position]);
    }
  }
  return inchikey_array;
}

function get_new_max_vals(data)
{
  var xy_min_max = [null,null,null,null];
  let x_values = [];
  let y_values = [];
  let x_col = [];
  let y_col = [];
  for(let i = 0; i < data.length;i++)
  {
    x_col = data[i].map(function(value,index){return value[0];});
    x_values.push.apply(x_values,x_col);
    y_col = data[i].map(function(value,index){return value[1];});
    y_values.push.apply(y_values,y_col);
  }
  xy_min_max[0] = Math.min.apply(null,x_values);
  xy_min_max[1] = Math.max.apply(null,x_values);
  xy_min_max[2] = Math.min.apply(null,y_values);
  xy_min_max[3] = Math.max.apply(null,y_values);

  x_dist = (xy_min_max[1] - xy_min_max[0])/5;
  y_dist = (xy_min_max[3] - xy_min_max[2])/2.5;

  xy_min_max[0] = Math.round(xy_min_max[0]-x_dist);
  xy_min_max[1] = Math.round(xy_min_max[1]+x_dist);
  xy_min_max[2] = Math.round(xy_min_max[2]-y_dist);
  xy_min_max[3] = Math.round(xy_min_max[3]+y_dist);

  return xy_min_max;
}

function extract_column(data,position)
{
  var separated_col = [];
  var concat_col = [];
  for(let i = 0; i < data.length; i++)
  {
    separated_col.push(null)
    separated_col[i]=data[i].map(function(value,index){return value[position];});
    concat_col.push.apply(concat_col,separated_col[i]);
  }
  return {separated_col:separated_col,concat_col:concat_col};
}

function get_scatter_option(data,chart_type,image,legend,tooltip_col)
{
  var get_rich_text = function(chart_data,image)
  {
    var rich_text = [];
    var curr_rich_text = {Attr:{width:40,align:'left'},Val:{width:60,align:'right'}};
    for(let i = 0; i < chart_data.length; i++)
    {
      if(data_size < data_size_threshold)
      {
        var image_url_collection = image[i];
        for(let j = 0; j < image_url_collection.length; j++)
        {
          var inchikey = chart_data[i][j][2].replace(/-/g,'');
          curr_rich_text[`${inchikey}`]={height:150,backgroundColor:{image:image_url_collection[j]}}
        }
      }
      rich_text.push(curr_rich_text);
      curr_rich_text = {Attr:{width:40,align:'left'},Val:{width:60,align:'right'}};
    }
    return rich_text;
  }
  var get_series = function(chart_data,chart_type,tooltip_col,rich_text)
  {
    var series_val = [];
    var choice_Min;
    var choice_diff;
    if(chart_type==2)
    {
      var choice_col = extract_column(chart_data,choice_param).concat_col;
      var choice_Max = Math.round(Math.max.apply(null,choice_col));
      choice_Min = Math.round(Math.min.apply(null,choice_col));
      choice_diff = choice_Max-choice_Min;   
    }
    for(let i = 0; i < chart_data.length; i++)
    {
      var series = {
        name: ChartOptions['legend'][i],
        data: chart_data[i],
        type: 'scatter',
        symbolSize: function(param){
          if(chart_type==1){return 15;}
          else{
          var minSize = 10;
          var maxSize = 40;
          var Sizediff = maxSize-minSize;
          var symbolSz = (param[choice_param]-choice_Min)*Sizediff/choice_diff+minSize;
          return symbolSz;}
        },
        emphasis: {
          focus: 'self',
          label: {
            show: true,
            formatter: function (param) {
              var return_list;
              var inchikey = param.data[2].replace(/-/g,'');
              if(data_size<data_size_threshold)
                {return_list = [`{${inchikey}|}`];}
              else
                {return_list = [];}
              for(let i = 0; i < tooltip_col.length; i++)
              {
                return_list.push(`{Attr|${tooltip_col[i]}}{Val|${param.data[i]}}`);
              }
              var rtext = return_list.join('\n');
              return rtext;
            },
            position: 'top',
            color: '#000',
            backgroundColor: '#fff',
            fontWeight: 'bold',
            padding: [3,4],
            width: 300,
            height: 100,
            fontSize: 10,
            align: 'center',
            verticalAlign: 'bottom',
            rich:rich_text[i],
          }
        },
      }
      series_val.push(series);
    }
    return series_val;
  }
  
  var rich_text = get_rich_text(data,image);
  var series_val = get_series(data,chart_type,tooltip_col,rich_text);
  var xy_min_max = get_new_max_vals(data);

  var option={
      legend: {
        right: '10%',
        top: '3%',
        data: legend,
      },
      xAxis: {
        type:'value',
        min: xy_min_max[0],
        max: xy_min_max[1],
        name: tooltip_col[0],
        nameLocation: 'start',
        nameTextStyle:{
          fontWeight:'bold',
          fontSize: 12,
          borderType: 'solid',},
        nameGap:30,},
      yAxis: {
        type:'value',
        min: xy_min_max[2],
        max: xy_min_max[3],
        name: tooltip_col[1],
        nameLocation: 'middle',
        nameTextStyle:{
          fontWeight:'bold',
          fontSize: 12,
          borderType: 'solid',},
        nameGap:50},
      series: series_val,
    };

  return option;
}

function process_histogram_data(data, num_of_bar)
{
  var extracted_col = extract_column(data,0);
  x_col = extracted_col.separated_col;
  x_values = extracted_col.concat_col;
  var x_min = Math.floor(Math.min.apply(null,x_values));
  var x_max = Math.floor(Math.max.apply(null,x_values))+1;
    
  let interval = (x_max-x_min)/num_of_bar;
  let threshold_arr = [];
  for(let i = 0; i <= num_of_bar; i++)
  {
    threshold_arr.push(x_min + i*interval);
  }
  let x_axis_val = [];
  for(let i = 0; i < num_of_bar; i++)
  {
    x_axis_val.push((threshold_arr[i]+threshold_arr[i+1])/2);
  }


  let count_arr = [];
  for(let i = 0; i < data.length; i++)
  {
    count_arr.push([]);
    for(let j = 0; j < num_of_bar; j++)
    {
      count_arr[i].push(0);
    }
  }
  for(let k = 0; k < data.length; k++)
  {
    for(let i = 0; i < x_col[k].length; i++)
    {
      for(let j = 0; j < num_of_bar; j++)
      {
        if(x_col[k][i] >= threshold_arr[j] && x_col[k][i] < threshold_arr[j+1])
       {
          count_arr[k][j]++;
        }
      }
    }
  }

  var aggregated_count = [];
  for(let i = 0; i< data.length; i++)
  {
    aggregated_count.push.apply(aggregated_count,count_arr[i]);
  }
  var y_max = Math.max.apply(null,aggregated_count);


  let processed_data = [];
  for(let i = 0; i < data.length; i++)
  {
    processed_data.push([]);
    for(let j = 0; j < num_of_bar; j++)
    {
      processed_data[i].push([x_axis_val[j],count_arr[i][j],Number(threshold_arr[j].toFixed(2)),Number((threshold_arr[j+1].toFixed(2)))]);
    }
  }

  return {processed_data:processed_data,x_min:x_min,x_max:x_max,y_max:y_max};
}



function get_histogram_option(data,legend,tooltip_col,num_of_bar)
{
  var histogram_data = process_histogram_data(data,num_of_bar);
  var series_val= []
  for(let i = 0; i < data.length; i++)
  {
    var series = {
      name: ChartOptions['legend'][i],
      data: histogram_data.processed_data[i],
      type: 'bar',
      barGap: '0%',
      emphasis:{
        focus: 'self',
        label: {
          show:true,
          color:'black',
          offset: [0,-100],
          fontSize: 15,
          backgroundColor: 'white',
          borderColor: 'black',
          borderWidth: 2,
          borderType: 'solid',
          padding:5,
          formatter: function(param){
            return [`Range [${param.data[2]},${param.data[3]}]`,`Count ${param.data[1]}`].join('\n');
          }
        }
      }
    }
    series_val.push(series);
  }

  option={
    legend: {
      right: '10%',
      top: '3%',
      data: ChartOptions['legend'],
    },
    xAxis: {
      type:'value',
      min: histogram_data.x_min,
      max: histogram_data.x_max,
      name: ChartOptions['header'][0],
      nameLocation: 'start',
      nameTextStyle:{
        fontWeight:'bold',
        fontSize: 12,
        borderType: 'solid',},
      nameGap:30,},
    yAxis: {
      type:'value',
      min: 0,
      max: histogram_data.y_max,
      name: 'Count',
      nameLocation: 'middle',
      nameTextStyle:{
        fontWeight:'bold',
        fontSize: 12,
        borderType: 'solid',},
      nameGap:50},
    series: series_val,
  };
  return option;
}

/*Echarts Lib, functional*/
var ChartOptions;
var Chart;
var option;
var chart_type;
var lexicographic_position;
var curr_num_of_bar = 20;
var legend;
var tooltip_col;
var curr_data;
var choice_param = 5;
var data_size;
var data_size_threshold = 1000;

var chartDom = document.getElementById('dynamic-chart');
if(chartDom) {
  Chart = echarts.init(chartDom,null);
  ChartOptions = JSON.parse(document.currentScript.nextElementSibling.textContent);
  tooltip_col = ChartOptions['header'];
  curr_data = ChartOptions['data'];
  var image = ChartOptions['image'];
  legend = ChartOptions['legend'];
  chart_type = ChartOptions['type'];
  lexicographic_position = ChartOptions['lexicographic_position'];
  data_size = ChartOptions['size'];
  
  if([1,2].includes(chart_type))
  {
    option = get_scatter_option(curr_data,chart_type,image,legend,tooltip_col);
  }
  else if([3].includes(chart_type))
  {
    option = get_histogram_option(curr_data,legend,tooltip_col,curr_num_of_bar);
  }
  option && Chart.setOption(option);

  var csv_input = document.getElementById('export-csv-val');
  csv_input.value = JSON.stringify(get_selected_inchikey(option,chart_type));
}

$(document).ready(function(){
  var num_filter = 0;
  var available_filter = []
  $('#add-filter').click(function(){
    var index = null;
    for(let i = 0; i <= num_filter;i++)
    {
      if(!available_filter.includes(i))
      {
        index = i;
        break;
      }
    }
    var new_div = document.createElement('div');
    new_div.id='filter_div_'+index;
    var new_select = document.createElement('select');
    new_select.id='filter_'+index;
    var header = ChartOptions['header'];
    var html = '';
    for(let i = 0;i<header.length;i++)
    {
      html += `<option value=${i}>${header[i]}</option>`;
    }
    var html_a="<a id='dynamic_input_" + index + "'>";
    html_a += ` ranges from <input type=text id=filter_input_${index}_0></input> to <input type=text id=filter_input_${index}_1></input>`;

    var button = document.createElement('button');
    button.id='remove_'+index;
    button.classList.add('filter-btn');
    const element = document.getElementById('filters');
    const elem_button = document.getElementById('add-filter');
    new_div.appendChild(button);
    new_div.appendChild(new_select);
    element.insertBefore(new_div,elem_button);
    $('#filter_'+index).append(html);
    $('#filter_div_'+index).append(html_a);
    num_filter++;
    available_filter.push(index);
    return;
  });
  $('#filters').on('change','select',function(events){
    var id = this.id;
    var split_id = id.split("_");
    var index = split_id[1];
    
    if(document.getElementById('dynamic_input_'+index))
    {
      $('#dynamic_input_'+index).remove();
    }
    var html="<a id='dynamic_input_" + index + "'>";
    if(this.value!=2 && this.value!=4)
    {
      html+=' ranges from ';
      html+=`<input type=text id=filter_input_${index}_0></input>`;
      html+=' to ';
      html+=`<input type=text id=filter_input_${index}_1></input>`;
    }
    else
    {
      html+=' contains ';
      html+=`<input type=text id=filter_input_${index}_0></input>`;
    }
    html+="</a>";
    $('#filter_div_'+index).append(html);
    return;
  })
  $('#filters').on('click','button.filter-btn',function(events){
    var id = this.id;
    var split_id = id.split("_");
    var deleteindex = split_id[1];
    $('#filter_div_'+deleteindex).remove();
    num_filter--;
    var new_arr = [];
    for(let i = 0;i<available_filter.length;i++)
    {
      if(available_filter[i]!=deleteindex)
      {
        new_arr.push(available_filter[i]);
      }
    }
    available_filter=new_arr;
    return;
  });
  $('#apply-filter').click(function(){
    curr_data = $.extend(true,{},ChartOptions['data']);
    if(available_filter.length==0)
    {
      console.log('Resetting');
    }
    else
    {
      var filter = null;
      var filter_val = [null,null];
      for(const slot of available_filter)
      {
        filter = document.getElementById('filter_'+slot);
        if(!lexicographic_position.includes[filter.value])
        {
          filter_val[0] = document.getElementById(`filter_input_${slot}_0`).value;
          filter_val[1] = document.getElementById(`filter_input_${slot}_1`).value;
          curr_data = numerical_filter(curr_data,filter.value,filter_val[0],filter_val[1]);
        }
        else
        {
          filter_val[0] = document.getElementById(`filter_input_${slot}_0`).value;
          curr_data = lexicographic_filter(curr_data,filter.value,filter_val[0]);
        }
      }
    }
    if(chart_type == 3)
    {
      var histogram_data = process_histogram_data(curr_data,curr_num_of_bar);
      for(const i of Object.keys(curr_data))
      {
        option['series'][i]['data'] = histogram_data.processed_data[i];
      }
      option['xAxis']['min']=histogram_data.x_min;
      option['xAxis']['max']=histogram_data.x_max;
      option['yAxis']['max']=histogram_data.y_max;
    }
    else if([1,2].includes(chart_type))
    {
      for(const i of Object.keys(curr_data))
      {
        option['series'][i]['data'] = curr_data[i];
      }
      var xy_min_max = get_new_max_vals(curr_data);
      option['xAxis']['min'] = xy_min_max[0];
      option['xAxis']['max'] = xy_min_max[1];
      option['yAxis']['min'] = xy_min_max[2];
      option['yAxis']['max'] = xy_min_max[3];
    }
    
    option && Chart.setOption(option);
    
    var csv_input = document.getElementById('export-csv-val');
    csv_input.value = JSON.stringify(get_selected_inchikey(option));
    return;
  });
  $('#apply-customization').click(function(){
    var all_options = document.getElementById('customize-options');
    var all_custom_select = all_options.querySelectorAll('select');
    for(var i = 0; i < all_custom_select.length; i++)
    {
      var select_field = all_custom_select[i];
      if(select_field.id == 'new-param-choice')
      {
        if(select_field.value != -1)
        {
          chart_type = 2;
          choice_param = select_field.value-1;
          option = get_scatter_option(curr_data,chart_type,image,legend,tooltip_col);
        }
        else
        {
          chart_type = 1;
          option = get_scatter_option(curr_data,chart_type,image,legend,tooltip_col);
        }
      }
    }
    var all_custom_input = all_options.querySelectorAll('input');
    for(var i = 0; i < all_custom_input.length; i++)
    {
      var input_field = all_custom_input[i];
      if(input_field.id == 'new-min-x' && checkValidFloat(input_field.value))
      {
        option['xAxis']['min'] = parseFloat(input_field.value);
      }
      if(input_field.id == 'new-max-x' && checkValidFloat(input_field.value))
      {
        option['xAxis']['max'] = parseFloat(input_field.value);
      }
      if(input_field.id == 'new-min-y' && checkValidFloat(input_field.value))
      {
        option['yAxis']['min'] = parseFloat(input_field.value);
      }
      if(input_field.id == 'new-max-y' && checkValidFloat(input_field.value))
      {
        option['yAxis']['max'] = parseFloat(input_field.value);
      }
      if(input_field.id == 'bin-num'  && checkValidInt(input_field.value))
      {
        curr_num_of_bar = parseInt(input_field.value);
        var histogram_data = process_histogram_data(curr_data,curr_num_of_bar);
        for(const i of Object.keys(curr_data))
        {
          option['series'][i]['data'] = histogram_data.processed_data[i];
        }
        option['xAxis']['min']=histogram_data.x_min;
        option['xAxis']['max']=histogram_data.x_max;
        option['yAxis']['max']=histogram_data.y_max;
      }
    }
    option && Chart.setOption(option);
  });
  $('#export-png').click(function(){
    var canvas_collection = document.getElementsByTagName('canvas');
    var dataURL = canvas_collection[0].toDataURL("image/png");
    var newTab = window.open('about:blank','image from canvas');
    newTab.document.write("<img src='"+dataURL+"' alt='from canvas'/>");
  });
  $('#export-sdf').click(function(){
    var canvas_collection = document.getElementsByTagName('canvas');
    var dataURL = canvas_collection[0].toDataURL("image/png");
    var newTab = window.open('about:blank','image from canvas');
    newTab.document.write("<img src='"+dataURL+"' alt='from canvas'/>");
  });
});

window.onresize = function() {
  Chart.resize();
};