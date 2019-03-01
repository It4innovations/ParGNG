#include "LinedList.h"


LinedList::LinedList(void)
{
	count = 0;
	head=NULL;
	tail=NULL;
	actual=NULL;
}


LinedList::~LinedList(void)
{
}
//Neuron* LinedList::FindBest(Neuron **neuronField)
//{
//	 double maxErrorSecond =-1;
//	 Neuron* temp;
//		for (data *i = tail; i->next!=NULL ; i=i->next)
//		{
//			if (neuronField[*(i->value)]->error > maxErrorSecond)
//                {
//                    maxErrorSecond = neuronField[*(i->value)]->error;
//                    temp = neuronField[*(i->value)];
//                }
//		}
//		return temp;
//}

void LinedList::ChangeItem(int * old,int * newItem, int newID)
{
	for (data *i = head; i!=NULL ; i=i->next)
	{
		if(i->value==old)
		{
			i->value=newItem;
			i->ID=newID;
			break;
		}
	}
}

void LinedList::IncreaseByOne(int limit)
{
		for (data *i = head; i!=NULL ; i=i->next)
		{
again:			if(*(i->value)==limit)
			{
				actual=i;
				i=i->next;
				DeleteActual();

				if(i==NULL)
					break;
				goto again;
			}else
			{
				*(i->value)+=1;
			}

		}
}

bool LinedList::Add(int * val, int IDval)
{
	data * temp=new data;
	temp->next=NULL;
	temp->value=val;
	temp->ID=IDval;

	if(head==NULL)
	{
		head=temp;
		tail=temp;
		temp->prev=NULL;
		
	}
	else
	{
		temp->prev=tail;
		tail->next=temp;
		tail=temp;
	}
	count++;
	return true;
}


bool LinedList::Delete(int * val)
{
	for (data *i = head; i!=NULL ; i=i->next)
	{
		if(i->value==val)
		{
			if((i->next==NULL) &&(i->prev==NULL))
			{
				head=NULL;
				tail=NULL;
			}
			else
			if(i->next==NULL)
			{
				i->prev->next=NULL;
				tail=i->prev;
			}else
				if(i->prev==NULL)
			{
				i->next->prev=NULL;
				head=i->next;
			}
			else
			{
								i->prev->next=i->next;
								i->next->prev=i->prev;
			}

			*(i->value)=-1;
			delete i;
			count--;
			break;
		}
	}
	return false;
}



bool LinedList::DeleteActual(void)
{
	if(actual!=NULL)
	{
			if((actual->next==NULL) &&(actual->prev==NULL))
			{
				head=NULL;
				tail=NULL;
			}
			else
			if(actual->next==NULL)
			{
				actual->prev->next=actual->next;
				tail=actual->prev;
			}else
				if(actual->prev==NULL)
			{
				actual->next->prev=actual->prev;
				head=actual->next;
			}
			else
			{
								actual->prev->next=actual->next;
								actual->next->prev=actual->prev;
			}
						*(actual->value)=-1;
			delete actual;
			count--;
	}
	actual=NULL;
	return false;
}

data* LinedList::Reset()
{
	actual=tail;
	return actual;
}

data* LinedList::Next()
{
	actual=actual->next;
	return actual;

}
int LinedList::size()
{
	return count;
}