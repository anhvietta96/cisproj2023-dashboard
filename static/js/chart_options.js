var set_data = JSON.parse(document.currentScript.nextElementSibling.textContent);

$(document).ready(function(){
    var num_set = 0;
    var available_set = []
    $('#add-set').click(function(){
      var available_slot = null;
      for(let i = 1; i <= num_set+1;i++)
      {
        if(!available_set.includes(i))
        {
          available_slot = i;
          break;
        }
      }

      var new_div = document.createElement('div');
      new_div.id='set_div_'+available_slot;
      var new_select = document.createElement('select');
      new_select.id='set_'+available_slot;
      new_select.name='set_'+available_slot;
      var html = '';

      for(let i = 0;i<set_data.length;i++)
      {
        html += `<option value=${i+1}>${set_data[i]}</option>`;
      }

      var button = document.createElement('button');
      button.id='remove_set_'+available_slot;
      button.classList.add('set-btn');
      button.type='button';

      const element = document.getElementById('set-selector');
      const elem_button = document.getElementById('add-set');
      new_div.appendChild(new_select);
      new_div.appendChild(button);
      element.appendChild(new_div);
      $('#set_'+available_slot).append(html);

      num_set++;
      available_set.push(available_slot);
      return;
    });

    $('#set-selector').on('click','button.set-btn',function(events){
      var id = this.id;
      var split_id = id.split("_");
      var deleteindex = split_id[2];
      $('#set_div_'+deleteindex).remove();
      num_set--;
      var new_arr = [];
      for(let i = 0;i<available_set.length;i++)
      {
        if(available_set[i]!=deleteindex)
        {
          new_arr.push(available_set[i]);
        }
      }
      available_set=new_arr;
      return;
    });
  });