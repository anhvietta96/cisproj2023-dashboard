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

function get_selected_inchikey(option)
{
  var inchikey_array = [];
  for(const group of Object.keys(option['series']))
  {
    console.log(group);
    for(let i = 0; i < option['series'][group]['data'].length; i++)
    {
      inchikey_array.push(option['series'][group]['data'][i][2]);
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
  let y_col = []
  for(let i = 0; i < data.length;i++)
  {
    x_col = data[i].map(function(value,index){return value[0];});
    x_values.push.apply(x_values,x_col);
    y_col = data[i].map(function(value,index){return value[1];});
    y_values.push.apply(y_values,y_col);
  }
  xy_min_max[0] = Math.round(Math.min.apply(null,x_values));
  xy_min_max[1] = Math.round(Math.max.apply(null,x_values));
  xy_min_max[2] = Math.round(Math.min.apply(null,y_values));
  xy_min_max[3] = Math.round(Math.max.apply(null,y_values));

  console.log(xy_min_max);

  x_dist = (xy_min_max[1] - xy_min_max[0])/5;
  y_dist = (xy_min_max[3] - xy_min_max[2])/2.5;

  xy_min_max[0] = xy_min_max[0]-x_dist;
  xy_min_max[1] = xy_min_max[1]+x_dist;
  xy_min_max[2] = xy_min_max[2]-y_dist;
  xy_min_max[3] = xy_min_max[3]+y_dist;

  return xy_min_max;
}

/*Echarts Lib, functional*/
var ChartOptions = null;
var Chart = null;
var option;

var chartDom = document.getElementById('dynamic-chart');
if(chartDom) {
  Chart = echarts.init(chartDom,null);
  ChartOptions = JSON.parse(document.currentScript.nextElementSibling.textContent);
  var tooltip_col = ChartOptions['header'];
  var chart_data = ChartOptions['data'];

  var xy_min_max = get_new_max_vals(chart_data);

  var series_val = [];
  var series = null;
  var rich_text = `{"Attr":{"width": 40,"align":"left"},"Val":{"width": 60,"align":"right"},`;
  var image_container = null;
  var inchikey = null;
  for(let i = 0;i<chart_data.length;i++)
  {
    image_container=[];
    image_url_collection = ChartOptions['image'][i];
    for(let j = 0;j<image_url_collection.length;j++)
    {
      image_container.push(null);
      image_container[j]=new Image();
      image_container[j].src=image_url_collection[j];
      inchikey = chart_data[i][j][2].replace(/-/g,'');
      rich_text += `"${inchikey}":{"height":150,"backgroundColor":{"image":"null"}},`
    }
    rich_text = rich_text.slice(0,-1) + '}';
    rich_text = JSON.parse(rich_text);
    for(let j = 0;j<image_url_collection.length;j++)
    {
      inchikey = chart_data[i][j][2].replace(/-/g,'');
      rich_text[`${inchikey}`]["backgroundColor"]["image"]=image_container[j];
    }
    

    series = {
      name: ChartOptions['name'][i],
      data: chart_data[i],
      type: 'scatter',
      symbolSize: 15,
      emphasis: {
        focus: 'self',
        label: {
          show: true,
          formatter: function (param) {
            inchikey = param.data[2].replace(/-/g,'');
            return_list = [`{${inchikey}|}`];
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
          rich:rich_text,
        }
      },
    }
    series_val.push(series);
    rich_text=`{"Attr":{"width": 60,"align":"left","padding":[0,10,0,10]},"Val":{"width": 90,"align":"right","padding":[0,10,0,10]},`
  }
  
  option={
    legend: {
      right: '10%',
      top: '3%',
      data: ChartOptions['legend'],
    },
    xAxis: {
      type:'value',
      min: xy_min_max[0],
      max: xy_min_max[1],
      name: ChartOptions['header'][0],
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
      name: ChartOptions['header'][1],
      nameLocation: 'middle',
      nameTextStyle:{
        fontWeight:'bold',
        fontSize: 12,
        borderType: 'solid',},
      nameGap:50},
    series: series_val,
  };
  option && Chart.setOption(option);

  var csv_input = document.getElementById('export-csv-val');
  csv_input.value = JSON.stringify(get_selected_inchikey(option));
  console.log(csv_input.value);
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
    var curr_data = $.extend(true,{},ChartOptions['data']);
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
        if(filter.value != 2 && filter.value != 4)
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
    for(const i of Object.keys(curr_data))
    {
      option['series'][i]['data'] = curr_data[i];
    }
    var xy_min_max = get_new_max_vals(curr_data);
    option['xAxis']['min'] = xy_min_max[0];
    option['xAxis']['max'] = xy_min_max[1];
    option['yAxis']['min'] = xy_min_max[2];
    option['yAxis']['max'] = xy_min_max[3];
    option && Chart.setOption(option);

    var csv_input = document.getElementById('export-csv-val');
    csv_input.value = JSON.stringify(get_selected_inchikey(option));
    return;
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